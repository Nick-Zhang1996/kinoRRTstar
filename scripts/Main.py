# test controller, fly polynomial trajectory from one point to another
import threading
import pickle
from threading import Event,Thread,Lock
from time  import time,sleep
from queue import Queue
from collections import namedtuple
from PidController import PidController
from scipy.spatial.transform import Rotation
from common import *
from math import degrees,cos,sin,radians
import matplotlib.pyplot as plt

import cflib
from cflib.crazyflie.log import LogConfig
#from cflib.crazyflie.syncLogger import SyncLogger
from cflib.crazyflie import Crazyflie

from Optitrack import Optitrack
from vicon import Vicon
from CommandTypes import *

from numeric_velocity import quad_fit_functional
from timeUtil import execution_timer

import os
import sys
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src/')
sys.path.append(base_dir)
from kinoRRT import *
from WorldVisualization import WorldVisualization
from cfcontroller import CFcontroller

from window import createWindow


class Main:
    def __init__(self,visual_tracker='vicon'):

        # if False block control from CF
        self.enable_control = True
        # common settings
        self.visual_tracker_freq = 120
        # crazyflie address
        self.uri = 'radio://0/80/2M/E7E7E7E7E7'

        # parameters to limit trajectory properties
        self.max_speed_limit = 3.0
        self.max_acc_limit = 10

        self.dt = 1.0/self.visual_tracker_freq
        # for simple performance tracking
        self.p = execution_timer(True)
        self.enable_log = Event()
        # flag for quitting
        self.quit_flag = Event()
        # keep track of spawned child threads
        self.child_threads = []

        # NOTE order of initialization may be important
        # Maybe put this in a separate function
        self.initKinoRrt()
        self.initExternalController()
        self.initLog()
        self.initVisualTracking(visual_tracker)
        self.initCrazyflie()
        self.initControllers()
        if (not self.enable_control):
            print_warning("control is disabled")

        return

    def initKinoRrt(self):
        # newly created window obstacle
        start_node = Node(-1.8, 0.6, -0.6)
        goal_node = Node(1.7,0,-2)

        world = World(-2,2,-0.7,1,-2.5,0)
        dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
        visual = WorldVisualization(dim)
        obstacles = []
        obstacles += createWindow((-1,0,0), 0.5, True, world)
        obstacles += createWindow((0.5,-1,0), 0.5, False, world)

        for obstacle in obstacles:
            world.addObstacle(obstacle)
            visual.addObstacle(obstacle)

        # TODO check start and goal are in world bobundary
        print_info("initializing kinoRRT*")
        self.rrt = rrt = KinoRrtStar(world, start_node, goal_node, 1800, 50)

        try:
            print_info("running kinoRRT*")
            t0 = time()
            rrt.run()
            elapsed = time()-t0
        except (RuntimeError):
            print_error("rrt error")

        print_info("rrt search finished in " + str(elapsed) +"sec")

        # sampled nodes
        critical_waypoint_n = rrt.prepareSolution()
        critical_waypoints = []
        for i in range(critical_waypoint_n):
            p = rrt.getNextWaypoint()
            assert (p.valid)
            critical_waypoints.append( (p.t,p.x,p.y,p.z) )
        # list of (t,x,y,z) 
        critical_waypoints = np.array(critical_waypoints)

        # waypoints, for accurate visualization
        self.traj_t = traj_t = rrt.getTrajectoryTime()
        tt = np.linspace(0,traj_t,1000)
        waypoints = []
        for this_t in tt:
            waypoint = rrt.getTrajectory(this_t)
            waypoints.append([waypoint.t, waypoint.x, waypoint.y, waypoint.z, waypoint.vx, waypoint.vy, waypoint.vz, waypoint.ax, waypoint.ay, waypoint.az])
        waypoints = np.array(waypoints)

        # print out critical information about trajectory
        max_speed = (np.max(waypoints[:,4]**2 + waypoints[:,5]**2 + waypoints[:,6]**2))**0.5
        max_acc = (np.max(waypoints[:,7]**2 + waypoints[:,8]**2 + waypoints[:,9]**2))**0.5
        print_info("total time : %.1f sec " %(traj_t))
        print_info("max speed : %.1f m/s " %(max_speed))
        print_info("max acc : %.1f m/s " %(max_acc))

        # rrt trajectory may demand unrealistic acceleration and velocity
        # rescale
        self.waypoints = waypoints
        self.processWaypoints()
        waypoints = self.waypoints


        ax = visual.visualizeWorld(show=False)
        ax.scatter(-critical_waypoints[:,1], critical_waypoints[:,2], -critical_waypoints[:,3], 'ro')
        ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'r')
        ax.set_xlabel("-x")
        ax.set_ylabel("y")
        ax.set_zlabel("-z")
        plt.show()

        print("start")
        print(waypoints[0])
        print("goal")
        print(waypoints[-1])

        response = input("press q + Enter to abort, Enter to execute")
        if (response == 'q'):
            print_warning("Aborting...")
            self.quit()
            exit(0)
            return
        print_info("saving trajectory")
        output = open("rrt_traj.p",'wb')
        pickle.dump(waypoints,output)
        output.close()
        print_ok("success")

        print_ok("Executing trajectory")

    # adjust time scale and spacial location
    def processWaypoints(self):
        waypoints = self.waypoints
        # scale time to match velocity to vehicle capabilities
        # find max speed and scale time respectively
        traj_t = waypoints[-1,0]
        max_speed = (np.max(waypoints[:,4]**2 + waypoints[:,5]**2 + waypoints[:,6]**2))**0.5
        max_acc = (np.max(waypoints[:,7]**2 + waypoints[:,8]**2 + waypoints[:,9]**2))**0.5

        speed_scale = self.max_speed_limit / max_speed
        acc_scale = self.max_acc_limit / max_acc
        self.time_scale = np.min([speed_scale, acc_scale,1.0])
        print_info("scaling factor = %.2f"%(self.time_scale))
        if (speed_scale < acc_scale):
            print_info("speed limited")
        else:
            print_info("acc limited")
        waypoints[:,0] /= self.time_scale
        waypoints[:,4:] *= self.time_scale

        max_speed = (np.max(waypoints[:,4]**2 + waypoints[:,5]**2 + waypoints[:,6]**2))**0.5
        max_acc = (np.max(waypoints[:,7]**2 + waypoints[:,8]**2 + waypoints[:,9]**2))**0.5

        print_info(" after scaling ")
        print_info("total time : %.1f sec " %(traj_t/self.time_scale))
        print_info("max speed : %.1f m/s " %(max_speed))
        print_info("max acc : %.1f m/s " %(max_acc))

        self.waypoints = waypoints

    def initExternalController(self):
        self.external_controller_t0 = None
        self.cfcontroller = CFcontroller(self.dt, self.time_scale, self.rrt)


    def initLog(self):
        # log vector: (x,y,z,rx,ry,rz,vx,vy,vz)
        self.log_vec = []
        # another format of log
        self.log_dict = {'thrust':0, 'target_vxy':(0,0)}
        return

    def initVisualTracking(self, visual_tracker):
        # local new state (from visual tracking) flag
        self.new_state = Event()
        self.visual_tracker = visual_tracker
        # state: (x,y,z,rx,ry,rz)
        self.drone_states = [0.0]*6
        self.drone_states_lock = Lock()
        # vel: (vx,vy,vz)
        self.drone_vel = [0.0,0.0,0.0]
        if (visual_tracker == 'optitrack'):
            self.vt = Optitrack(freq = self.visual_tracker_freq)
            # optitrack Id for drone
            # this is the id listed in Motive software
            self.vt_id = 6
            # Optitrack interface has a dynamically assigned internal Id that differ from global optitrack objet ID
            # the internal id is used to retrieve state from optitrack instance
            self.optitrack_internal_id = self.vt.getInternalId(self.optitrack_id)
        elif (visual_tracker == 'vicon'):
            #vicon item id
            self.vt = Vicon(daemon=True)
            self.vt.getViconUpdate()
            sleep(0.1)
            self.vt_id = self.vt.getItemID('nick_cf')

        # start daemon thread
        if (self.visual_tracker == 'optitrack'):
            self.child_threads.append(Thread(target=self.optitrackUpdateThread))
            self.child_threads[-1].start()
            print_info("starting thread optitrackUpdateThread")
        elif (self.visual_tracker == 'vicon'):
            self.child_threads.append(Thread(target=self.viconUpdateThread))
            self.child_threads[-1].start()
            print_info("starting thread viconUpdateThread")
        return

    def initControllers(self):
        dt = self.dt
        # indicate whether a new command is available
        self.new_command = Event()
        # commands pipeline
        # Vel or Pos or Planar
        self.commands = Queue()

        self.x_pids = PidController(2,0,0,dt,0,20)
        self.y_pids = PidController(2,0,0,dt,0,20)
        self.z_pids = PidController(2,0,0,dt,0,20)
        self.vx_pids = PidController(25,1,1,dt,5,30)
        self.vy_pids = PidController(25,1,1,dt,5,30)
        self.vz_pids = PidController(25,15,1,dt,1,30)
        self.yaw_pid = PidController(2,0,0,dt,0,20)

        # output satuation for position controller
        self.vxy_limit = 1.0
        self.vz_limit = 1.0

        self.baseThrust = 42500
        self.minThrust = 20000
        self.thrustScale = 1000.0

        # start PID controller thread
        self.external_controller_active = Event()
        self.child_threads.append(Thread(target=self.pidCrazyflieControl))
        self.child_threads[-1].start()
        print_info("starting thread crazyflieControl")

        self.child_threads.append(Thread(target=self.externalCrazyflieControl))
        self.child_threads[-1].start()
        print_info("starting thread ExternalCrazyflieControl")

    def initCrazyflie(self):
        cflib.crtp.init_drivers(enable_debug_driver=False)
        uri = self.uri
        self.connect_signal = Event()
        self.cf = Crazyflie(rw_cache='.cache')
        self.cf.connected.add_callback(self.connectedCallback)
        self.cf.disconnected.add_callback(self.disconnectedCallback)
        self.cf.open_link(uri)
        self.connect_signal.wait()

    def issueCommand(self, command):
        self.commands.put(command)
        self.new_command.set()
        return

    # main entry point
    def run(self,):

        # guide crazyflie to initial position
        print_ok("taking off")
        (x,y,z,_,_,_) = self.drone_states
        print_info("current pos")
        print_info(x,y,z)
        self.issueCommand(Planar(0,0,-0.3))
        response = input("press Enter to continue(go to start pos), q+enter to quit \n")
        if (response == 'q'):
            print_warning("Aborting...")
            self.quit()
            exit(0)
            return

        print_ok("going to start position")
        retval = self.cfcontroller.getTrajectory(0.0)
        target_x,target_y,target_z = retval
        print_info("start pos")
        print_info(retval)
        cmd = Pos(x=target_x, y=target_y, z=target_z)
        self.issueCommand(cmd)
        response = input("press Enter to continue, q+enter to quit \n")
        if (response == 'q'):
            print_warning("Aborting...")
            self.quit()
            exit(0)
            return

        print("External Control Starts")
        self.enable_log.set()
        self.external_controller_t0 = time()
        self.external_controller_active.set()

        input("press Enter to land (and stop log) \n")
        print_ok("landing")
        self.issueCommand(Planar(0,0,-0.1))
        self.external_controller_active.clear()
        self.enable_log.clear()

        input("press Enter to shutdown \n")
        self.quit()
        return

    def quit(self,):
        self.quit_flag.set()
        print_info("quit flag set")

        print_info("joining child threads...")
        for p in self.child_threads:
            p.join()
        print_info("all joined ")

        self.vt.quit()
        self.logFilename = "./log.p"
        output = open(self.logFilename,'wb')
        pickle.dump(self.log_vec,output)
        output.close()
        self.cf.close_link()

    def optitrackUpdateThread(self):
        # update state
        while not self.quit_flag.isSet():
            # wait for optitrack to get new state update
            ret = self.vt.newState.wait(0.1)
            if (not ret):
                continue
            lock = self.drone_states_lock.acquire(timeout=0.01)
            if lock:
                self.vt.newState.clear()
                (x,y,z,rx,ry,rz) = state = self.vt.getLocalState(self.optitrack_internal_id)
                # TODO verify
                if (self.vt.lost[0].isSet()):
                    self.quit_flag.set()
                    print_warning("Lost track of object")

                self.drone_states = state
                self.drone_vel = self.vt.getLocalVelocity(self.optitrack_internal_id)
                #self.log_vec.append(self.drone_states+tuple(self.drone_vel)+(self.log_dict['thrust'],self.log_dict['target_vxy']))
                if (self.enable_log.is_set()):
                    log_entry = (time(),) + tuple(self.drone_states)
                    self.log_vec.append(log_entry)
                self.new_state.set()
                self.drone_states_lock.release()

    def viconUpdateThread(self):
        # update state
        while not self.quit_flag.isSet():
            ret = self.vt.newState.wait(0.1)
            if (not ret):
                continue
            lock = self.drone_states_lock.acquire(timeout=0.01)
            if lock:
                self.vt.newState.clear()
                self.drone_states = (x,y,z,rx,ry,rz) = state = self.vt.getState(self.vt_id)
                self.drone_vel = (vx,vy,vz) = self.vt.getVelocity(self.vt_id)
                #print("%7.3f, %7.3f, %7.3f " %(vx,vy,vz))
                if (self.enable_log.is_set()):
                    log_entry = (time(),) + tuple(self.drone_states)
                    self.log_vec.append(log_entry)
                self.new_state.set()
                self.drone_states_lock.release()

    def connectedCallback(self, uri):
        self.connect_signal.set()
        print_ok("Connection to Crazyflie established.")

    def disconnectedCallback(self, uri):
        print_warning("Crazyflie disconnected")

    def activate_high_level_commander(self,scf):
        try:
            scf.cf.param.set_value('commander.enHighLevel', '1')
        except Exception as e:
            print_error(e)
            return

    # read new command and call requestXXX() functions
    def externalCrazyflieControl(self):
        try:
            cf = self.cf

            while not self.quit_flag.isSet():
                if (not self.external_controller_active.is_set() ):
                    sleep(0.01)
                    continue
                (x,y,z,rx,ry,rz) = self.drone_states
                (vx,vy,vz) = self.drone_vel

                if (z<-4.0):
                    print_warning(" exceeding maximum allowable height ")
                    print_info("switching to safety mode")
                    self.issueCommand(Planar(0,0,-0.1))
                    self.external_controller_active.clear()
                    return

                drone_state = (x,y,z,vx,vy,vz,rx,ry,rz)
                ret = self.cfcontroller.control(time()-self.external_controller_t0, drone_state)
                if (ret is None):
                    print_info("external control finished, yielding control")
                    self.issueCommand(Planar(0,0,-0.1))
                    self.external_controller_active.clear()
                    self.enable_log.clear()
                else:
                    (target_roll_deg, target_pitch_deg, target_yawrate_deg_s, target_thrust_raw) = ret

                    target_thrust_raw = int(np.clip(target_thrust_raw,self.minThrust,0xFFFF))

                    if (self.enable_control):
                        cf.commander.send_setpoint(target_roll_deg,-target_pitch_deg,target_yawrate_deg_s,target_thrust_raw)

        except Exception as e:
            print_error("crazyflieControl: "+str(e))
        return

    # read new command and call requestXXX() functions
    def pidCrazyflieControl(self):
        try:
            cf = self.cf
            # an empty command must be sent to unlock thrust protection
            cf.commander.send_setpoint(0,0,0,0)
            self.log_dict['thrust'] = 0
            command = None

            while not self.quit_flag.isSet():
                if ( self.external_controller_active.is_set() ):
                    sleep(0.01)
                    continue
                # TODO add failsafe
                s = self.drone_states
                #print("%.2f, %.2f %.2f ::: %.2f, %.2f %.2f"%(s[0],s[1],s[2],degrees(s[3]),degrees(s[4]),degrees(s[5])))
                vel = self.drone_vel
                #print("%.2f, %.2f %.2f "%(vel[0],vel[1],vel[2]))

                ret = self.new_command.wait(0.02)
                if ret:
                    #print("New command received")
                    self.new_command.clear()
                    candidate_command = self.commands.get()
                    if candidate_command is None:
                        print_warning("expecting new command, received None")
                    else:
                        command = candidate_command

                if isinstance(command,Vel):
                    #print("Velocity Command")
                    try:
                        self.requestVelocity(cf,(command.vx, command.vy, command.vz))
                    except Exception as e:
                        print_error("requestVelocity: "+str(e))
                elif isinstance(command,Pos):
                    #print("Position Command")
                    try:
                        self.requestPosition(cf,(command.x, command.y, command.z))
                    except Exception as e:
                        print_error("requestPosition: "+str(e))
                elif isinstance(command,Planar):
                    #print("Planar Command")
                    try:
                        self.requestPlanar(cf,(command.vx, command.vy, command.z))
                    except Exception as e:
                        print_error("requestPlanar: "+str(e))
                else:
                    pass
                    #print("unknown command")

            # stop everything before returning
            print_info("sending zero thrust command to CF")
            cf.commander.send_setpoint(0,0,0,0)
            self.log_dict['thrust'] = 0
            sleep(0.1)
        except Exception as e:
            print_error("crazyflieControl: "+str(e))
        return

    def requestPlanar(self,cf,planar_command):
        target_vx,target_vy,target_z = planar_command
        target_vz = self.z_pids.control(planar_command[2],self.drone_states[2])
        self.requestVelocity(cf, (target_vx, target_vy, target_vz))
        return


    def requestPosition(self,cf,pos_command):
        target_x,target_y,target_z = pos_command
        target_vx_world = self.x_pids.control(target_x, self.drone_states[0])
        target_vy_world = self.y_pids.control(target_y, self.drone_states[1])
        target_vz_world = self.z_pids.control(target_z, self.drone_states[2])

        target_vx_world = np.clip(target_vx_world, -self.vxy_limit,self.vxy_limit)
        target_vy_world = np.clip(target_vy_world, -self.vxy_limit,self.vxy_limit)
        target_vz_world = np.clip(target_vz_world, -self.vz_limit,self.vz_limit)

        target_v = (target_vx_world,target_vy_world,target_vz_world)
        self.requestVelocity(cf,target_v)
        return

    def requestVelocity(self,cf,vel_command, yaw_des = 0.0):
        # convert velocity to vehicle frame
        (x,y,z,rx,ry,rz) = state = self.drone_states
        #print(x,y,z,degrees(rz))
        # TODO optimize
        r = Rotation.from_euler("Z",[rz],degrees=False)

        try:
            target_v_local = r.inv().apply(vel_command).flatten()
            actual_v_local = r.inv().apply(np.array(self.drone_vel).flatten()).flatten()
            #print("%.2f, %.2f, %.2f"%(degrees(rx),degrees(ry),degrees(rz)))
            #print("%.2f, %.2f, %.2f"%(actual_v_local[0],actual_v_local[1],actual_v_local[2]))
            #print("%.2f, %.2f, %.2f"%(self.drone_vel[0],self.drone_vel[1],self.drone_vel[2]))

            # in crazyflie's ref frame roll to right is positive
            # in deg
            target_pitch_deg = -self.vx_pids.control(target_v_local[0], actual_v_local[0])
            target_pitch_deg = np.clip(target_pitch_deg, -30, 30)

            target_roll_deg = self.vy_pids.control(target_v_local[1], actual_v_local[1])
            target_roll_deg = np.clip(target_roll_deg, -30, 30)

            # NOTE assume z to point downward, NED frame
            #print("target: %.2f, actual %.2f"%(target_v_local[2],self.drone_vel[2]))
            # NOTE using world frame z velocity
            target_thrust = self.baseThrust - self.vz_pids.control(target_v_local[2], self.drone_vel[2]) * self.thrustScale
            target_thrust = int(np.clip(target_thrust,self.minThrust,0xFFFF))

            yaw_diff = yaw_des - rz
            yaw_diff = (yaw_diff + np.pi)%(2*np.pi) - np.pi
            target_yawrate_deg_s = self.yaw_pid.control(degrees(rz + yaw_diff), degrees(rz))



            #print_warning(" crazyflie command blocked ")
            if (self.enable_control):
                cf.commander.send_setpoint(target_roll_deg,-target_pitch_deg,target_yawrate_deg_s,target_thrust)
            #print_info("sending command %.2f %.2f %.2f %d"%(target_roll_deg,-target_pitch_deg,-target_yawrate_deg_s,target_thrust))
            self.log_dict['thrust'] = target_thrust
            #print(target_roll_deg,-target_pitch_deg,target_thrust)
            #print("thrust = %d"%target_thrust)
        except Exception as e:
            print_error(e)
        #print("\t roll: %.1f, pitch: %.1f, yawrate: %.1f, thrust %d"%(target_roll_deg,target_pitch_deg,target_yawrate_deg_s,target_thrust))
        return

    # unused
    def wait_for_param_download(self,scf):
        while not scf.cf.param.is_updated:
            time.sleep(1.0)
        print('Parameters downloaded for', scf.cf.link_uri)

    # return a tuple of (x,y,z), parameterized by time
    def getWaypointByTime(self,t):
        tf = self.waypoints[-1,0]
        if (t < tf):
            # find closest waypoint
            #i = np.argmin(np.abs(t - self.waypoints[:,0]))

            # note: t_i-1 < t <= t_i
            i = np.searchsorted( self.waypoints[:,0],t, side='right')
            # interpolate between two points
            alfa = (t - self.waypoints[i-1,0]) - (self.waypoints[i,0] - self.waypoints[i-1,0])
            interpolated = (1-alfa) * self.waypoints[i-1,1:] + alfa * self.waypoints[i,1:]
            return interpolated
        else:
            return None

    # get current vehicle state
    # send to velocity planner
    # call velocity controller
    # update velocity commands
    def controlThread(self,):
        t0 = time()

        while not self.quit_flag.isSet():
            ret = self.new_state.wait(0.1)
            if (not ret):
                continue
            self.new_state.clear()

            retval = self.getWaypointByTime(time()-t0)
            if (retval is None):
                print_info("control sequence complete")
                # TODO add elegant start
                cmd = Planar(0,0,-0.1)
                self.commands.put(cmd)
                self.new_command.set()
                return
            target_x,target_y,target_z = retval
            cmd = Pos(x=target_x, y=target_y, z=target_z)
            self.commands.put(cmd)
            self.new_command.set()
            sleep(0.05)



if __name__ == '__main__':
    ins = Main()
    ins.run()
    print_info("program finish")
