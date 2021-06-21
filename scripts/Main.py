# fly multiple drone
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

import cflib
from cflib.crazyflie.log import LogConfig
#from cflib.crazyflie.syncLogger import SyncLogger
from cflib.crazyflie import Crazyflie

from Optitrack import Optitrack
from vicon import Vicon
from CommandTypes import *

from numeric_velocity import quad_fit_functional
from PlanarController import planarController,planarControllerExample
from timeUtil import execution_timer

import os
import sys
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src/')
sys.path.append(base_dir)
from kinoRRT import *
from WorldVisualization import WorldVisualization


class Main:
    def __init__(self,visual_tracker='vicon'):

        # common settings put here for easy access
        self.visual_tracker_freq = 120
        self.dt = 1.0/self.visual_tracker_freq
        # crazyflie address
        self.uri = 'radio://0/80/2M/E7E7E7E7E7'

        # for simple performance tracking
        self.p = execution_timer(True)
        # flag for quitting
        self.quit_flag = Event()
        # keep track of spawned child threads
        self.child_threads = []

        # NOTE order of initialization may be important
        # Maybe put this in a separate function
        self.initLog()
        self.initVisualTracking(visual_tracker)
        self.initControllers()
        self.initCrazyflie()
        self.initKinoRrt()

        return

    def initKinoRrt(self):
        world = World(20,10,10)
        visual = WorldVisualization(size=(world.x,world.y,world.z))

        obs1 = Box(6,0,0,6+2,5,10)
        obs2 = Box(6,5,0,6+2,5+5,5)
        obs3 = Box(12,5,0,12+2,5+5,10)
        obs4 = Box(12,0,5,12+2,5,5+5)

        obstacles = [obs1, obs2, obs3, obs4]
        for obstacle in obstacles:
            world.addObstacle(obstacle)
            visual.addObstacle(obstacle)


        start_node = Node(2, 2, 2)
        goal_node = Node(18, 8, 8)

        rrt = KinoRrtStar(world, start_node, goal_node, 600, 10)

        try:
            t0 = time()
            rrt.run()
            elapsed = time()-t0
        except (RuntimeError):
            print("rrt error")
        print("rrt search finished in " + str(elapsed) +"sec")

        waypoint_n = rrt.prepareSolution()

        waypoints = []
        for i in range(waypoint_n):
            p = rrt.getNextWaypoint()
            assert (p.valid)
            waypoints.append( (p.x,p.y,p.z) )

        # list of (x,y,z) in reverse order
        waypoints = np.array(waypoints[::-1])
        #print("waypoints")
        #print(waypoints)

        ax = visual.visualizeWorld(show=False)
        #ax.plot(waypoints[:,0], waypoints[:,1], waypoints[:,2],'b')
        ax.scatter(waypoints[:,0], waypoints[:,1], waypoints[:,2], 'ro')
        plt.show()

        response = input("press q + Enter to abort, Enter to execute")
        if (response == 'q'):
            print_warning("Aborting...")
            self.quit()
            exit(0)
            return
        print_ok("Executing trajectory")

        # calculate a reasonable speed
        self.waypoints = waypoints


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

        # output satuation for position controller
        self.vxy_limit = 1.0
        self.vz_limit = 1.0

        self.baseThrust = 42500
        self.minThrust = 20000
        self.thrustScale = 1000.0

        # start controller thread
        self.child_threads.append(Thread(target=self.crazyflieControl))
        self.child_threads[-1].start()
        print_info("starting thread crazyflieControl")

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
        print_info("taking off")
        self.issueCommand(Planar(0,0,-0.3))
        sleep(1.0)

        print_info("going to start position")
        retval = self.getWaypointByTime(0.0)
        target_x,target_y,target_z = retval
        cmd = Pos(x=target_x, y=target_y, z=target_z)
        self.issueCommand(cmd)
        sleep(3.0)

        print("Main Control Starts")
        p = threading.Thread(target=self.controlThread)
        self.child_threads.append(p)
        self.child_threads[-1].start()

        input("press Enter to stop")
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
                log_entry = (time(),) + tuple(np.drone_states)
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
                self.log_vec.append(self.drone_states)
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
    def crazyflieControl(self):
        try:
            cf = self.cf
            # an empty command must be sent to unlock thrust protection
            cf.commander.send_setpoint(0,0,0,0)
            self.log_dict['thrust'] = 0
            command = None

            while not self.quit_flag.isSet():
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

    def requestVelocity(self,cf,vel_command):
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
            target_yawrate_deg_s = 0
            #print_warning(" crazyflie command blocked ")
            cf.commander.send_setpoint(target_roll_deg,-target_pitch_deg,-target_yawrate_deg_s,target_thrust)
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
        # parameters for a circle
        T = 10.0
        tf = 10.0
        # trajectory: sin(2*np.pi/T * t), cos(2*np.pi/T * t), -0.3
        # start at t=0 -> (0, 1, -0.3)
        if (t < tf):
            return (sin(2*np.pi/T * t), cos(2*np.pi/T * t), -0.3)
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
                sleep(1)
                self.cf.commander.send_setpoint(0,0,0,0)
                sleep(0.01)
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
