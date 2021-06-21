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

import cflib.crtp
from cflib.crazyflie.log import LogConfig
from cflib.crazyflie.swarm import CachedCfFactory
from cflib.crazyflie.swarm import Swarm
from cflib.crazyflie.syncLogger import SyncLogger

from Optitrack import Optitrack

from numeric_velocity import quad_fit_functional
from PlanarController import planarController,planarControllerExample

# command types
class Vel:
    def __init__(self,vx=0,vy=0,vz=0):
        self.vx = vx
        self.vy = vy
        self.vz = vz
        return

class Pos:
    def __init__(self,x=0,y=0,z=0):
        self.x = x
        self.y = y
        self.z = z
        return

class Planar:
    def __init__(self,vx=0,vy=0,z=0):
        self.vx = vx
        self.vy = vy
        self.z = z
        return

class MultipleDrone:
    def __init__(self,):
        self.optitrack_freq = 50

        # drone 1,2,3 ------
        self.uris = [
            'radio://0/80/2M/E7E7E7E701',  # cf_id 0
            'radio://0/80/2M/E7E7E7E702',  # cf_id 1
            'radio://0/80/2M/E7E7E7E703',  # cf_id 3
        ]

        # optitrack Id for each drone, id0, id1, id2
        # this is the id listed in Motive software
        self.optitrack_ids = [3,4,5]

        # drone 1,2
        self.uris = [
            'radio://0/80/2M/E7E7E7E701',  # cf_id 0
            'radio://0/80/2M/E7E7E7E702',  # cf_id 1
        ]

        # optitrack Id for each drone, id0, id1, id2
        # this is the id listed in Motive software
        self.optitrack_ids = [3,4,5]

        # target altitude
        #self.target_z = -0.3
        # target altitude during maneuver for each drone
        self.target_z_array = [-0.3, -0.3, -0.3]

        # set planar controller
        self.multiAgentControl = planarController

        # local new state indicator
        self.new_state = Event()

        # flag for quitting
        self.quit_flag = Event()
        # spawned child threads
        self.child_threads = []
        self.drone_count = len(self.uris)
        # commands
        # Vel or Pos or Planar
        self.commands = [Queue() for _ in range(len(self.uris))]

        # state: (x,y,z,rx,ry,rz)
        self.drone_states = [[0]*6 for _ in range(len(self.uris))]
        self.drone_states_lock = Lock()
        # vel: (vx,vy,vz)
        self.drone_vel = [[0,0,0] for _ in range(len(self.uris))]
        # indicate whether a new command is available
        self.new_command = [Event() for _ in range(len(self.uris))]
        # for each drone, (x,y,z,rx,ry,rz,vx,vy,vz)
        self.log_vec = [[] for _ in range(len(self.uris))]
        self.log_dict = {'thrust':0, 'target_vxy':(0,0)}

        dt = 1.0/self.optitrack_freq
        self.x_pids = [PidController(2,0,0,dt,0,20) for _ in range(len(self.uris))]
        self.y_pids = [PidController(2,0,0,dt,0,20) for _ in range(len(self.uris))]
        self.z_pids = [PidController(2,0,0,dt,0,20) for _ in range(len(self.uris))]
        self.vx_pids = [PidController(25,1,1,dt,5,30) for _ in range(len(self.uris))]
        self.vy_pids = [PidController(25,1,1,dt,5,30) for _ in range(len(self.uris))]
        self.vz_pids = [PidController(25,15,1,dt,1,30) for _ in range(len(self.uris))]

        self.vxy_limit = 1.0
        self.vz_limit = 1.0

        self.baseThrust = 42500
        self.minThrust = 20000
        self.thrustScale = 1000.0
        self.vel_est_n_step = 5
        self.xyz_history = []
        self.quad_fit = quad_fit_functional(self.vel_est_n_step, dt)

        self.op = Optitrack(freq = self.optitrack_freq)
        # Optitrack interface has a dynamically assigned internal Id that differ from global optitrack objet ID
        # the internal id is used to retrieve state from optitrack instance
        self.optitrack_internal_ids = [self.op.getInternalId(i) for i in self.optitrack_ids]
        return

    def quit(self,):
        self.quit_flag.set()
        for p in self.child_threads:
            p.join()
        self.op.quit()
        self.logFilename = "./log.p"
        output = open(self.logFilename,'wb')
        pickle.dump(self.log_vec,output)
        output.close()

    def optitrackUpdateThread(self):
        # update state
        while not self.quit_flag.isSet():
            # wait for optitrack to get new state update
            ret = self.op.newState.wait(0.1)
            if (not ret):
                continue
            self.op.newState.clear()
            with self.drone_states_lock:
                for i in range(self.drone_count):
                    (x,y,z,rx,ry,rz) = state = self.op.getLocalState(self.optitrack_internal_ids[i])
                    if (self.op.lost[i].isSet()):
                        self.quit_flag.set()
                        print_warning("Lost track of object")

                    self.drone_states[i] = state
                    self.drone_vel[i] = self.op.getLocalVelocity(self.optitrack_internal_ids[i])
                    self.log_vec[i].append(self.drone_states[i]+tuple(self.drone_vel[i])+(self.log_dict['thrust'],self.log_dict['target_vxy'][0],self.log_dict['target_vxy'][1],))
            self.new_state.set()

    def swarmThread(self,swarm):
        swarm.parallel_safe(self.crazyflieControl)

    def run(self,):
        self.child_threads.append(Thread(target=self.optitrackUpdateThread))
        self.child_threads[-1].start()

        cflib.crtp.init_drivers(enable_debug_driver=False)
        factory = CachedCfFactory(rw_cache='./cache')
        with Swarm(self.uris, factory=factory) as swarm:

            swarm.parallel(self.wait_for_param_download)
            #swarm.parallel_safe(self.activate_high_level_commander)

            print("starting swarm...")
            thread = Thread(target=self.swarmThread,args=(swarm,))
            thread.start()

            # guide each crazyflie to initial position
            # go to different height
            print("taking off")
            self.commands[0].put(Planar(0,0,-0.3))
            self.commands[1].put(Planar(0,0,-0.3))
            #self.commands[2].put(Planar(0,0,-0.3))
            self.new_command[0].set()
            self.new_command[1].set()
            #self.new_command[2].set()
            sleep(3)

            print("Going to initial positions")
            # for circle
            #self.commands[0].put(Pos(0.6,0,-0.3))
            #self.commands[1].put(Pos(0.3,0,-0.3))
            # for multiagent controller, 3 agents
            '''
            self.commands[0].put(Pos(1/3.0,1/15.0,-0.3))
            self.commands[1].put(Pos(2/3.0,1/3.0,-0.4))
            self.commands[2].put(Pos(1/3.0,4/15.0,-0.5))
            self.new_command[0].set()
            self.new_command[1].set()
            self.new_command[2].set()
            '''
            # for multiagent controller, 2 agents
            self.commands[0].put(Pos(0.55/1.5,0.7/1.5,-0.3))
            self.commands[1].put(Pos(-0.25/1.5,-0.3/1.5,-0.3))
            self.new_command[0].set()
            self.new_command[1].set()
            sleep(3)

            '''
            self.commands[0].put(Planar(0,0,-0.0))
            self.commands[1].put(Planar(0,0,-0.0))
            self.commands[2].put(Planar(0,0,-0.0))
            self.new_command[0].set()
            self.new_command[1].set()
            self.new_command[2].set()

            '''
            print("Main Control Starts")
            p = threading.Thread(target=self.controlThread)
            self.child_threads.append(p)
            self.child_threads[-1].start()
            input("press Enter to stop")

            self.quit()
        return

    def activate_high_level_commander(self,scf):
        try:
            scf.cf.param.set_value('commander.enHighLevel', '1')
        except Exception as e:
            print_error(e)
            return

    # execute the velocity commands for each crazyflie
    def crazyflieControl(self,scf):
        try:
            cf = scf.cf
            this_id = self.uris.index(cf.link_uri)
            #print("control thread started id:%d"%this_id)
            # an empty command must be sent to unlock thrust protection
            cf.commander.send_setpoint(0,0,0,0)
            self.log_dict['thrust'] = 0
            command = None

            while not self.quit_flag.isSet():
                # TODO add failsafe
                s = self.drone_states[this_id]
                #print("%.2f, %.2f %.2f ::: %.2f, %.2f %.2f"%(s[0],s[1],s[2],degrees(s[3]),degrees(s[4]),degrees(s[5])))
                vel = self.drone_vel[this_id]
                #print("%.2f, %.2f %.2f "%(vel[0],vel[1],vel[2]))

                ret = self.new_command[this_id].wait(0.02)
                if ret:
                    #print("id:%d New command received"%this_id)
                    self.new_command[this_id].clear()
                    candidate_command = self.commands[this_id].get()
                    if candidate_command is None:
                        print_warning("received none as command")
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
            print("id: %d stopping..."%(this_id))
            cf.commander.send_setpoint(0,0,0,0)
            self.log_dict['thrust'] = 0
            sleep(0.1)
        except Exception as e:
            print_error("crazyflieControl: "+str(e))
        return

    def requestPlanar(self,cf,planar_command):
        this_id = self.uris.index(cf.link_uri)
        target_vx,target_vy,target_z = planar_command
        target_vz = self.z_pids[this_id].control(planar_command[2],self.drone_states[this_id][2])
        self.requestVelocity(cf, (target_vx, target_vy, target_vz))
        return


    def requestPosition(self,cf,pos_command):
        this_id = self.uris.index(cf.link_uri)

        target_x,target_y,target_z = pos_command
        target_vx_world = self.x_pids[this_id].control(target_x,self.drone_states[this_id][0])
        target_vy_world = self.y_pids[this_id].control(target_y,self.drone_states[this_id][1])
        target_vz_world = self.z_pids[this_id].control(target_z,self.drone_states[this_id][2])

        target_vx_world = np.clip(target_vx_world,-self.vxy_limit,self.vxy_limit)
        target_vy_world = np.clip(target_vy_world,-self.vxy_limit,self.vxy_limit)
        target_vz_world = np.clip(target_vz_world,-self.vz_limit,self.vz_limit)

        target_v = (target_vx_world,target_vy_world,target_vz_world)
        self.requestVelocity(cf,target_v)
        return

    def requestVelocity(self,cf,vel_command):
        this_id = self.uris.index(cf.link_uri)
        # convert velocity to vehicle frame
        (x,y,z,rx,ry,rz) = state = self.drone_states[this_id]
        r = Rotation.from_euler("Z",[rz],degrees=False)

        try:
            target_v_local = r.inv().apply(vel_command).flatten()
            actual_v_local = r.inv().apply(np.array(self.drone_vel[this_id]).flatten()).flatten()
            #print("%.2f"%(degrees(rz)))
           # print("%.2f, %.2f, %.2f"%(degrees(rx),degrees(ry),degrees(rz)))


            #print("%.2f, %.2f, %.2f"%(actual_v_local[0],actual_v_local[1],actual_v_local[2]))
            #print("%.2f, %.2f, %.2f"%(self.drone_vel[this_id][0],self.drone_vel[this_id][1],self.drone_vel[this_id][2]))


            # in crazyflie's ref frame roll to right is positive
            # in deg
            target_pitch_deg = -self.vx_pids[this_id].control(target_v_local[0], actual_v_local[0])
            target_pitch_deg = np.clip(target_pitch_deg, -30, 30)

            target_roll_deg = self.vy_pids[this_id].control(target_v_local[1], actual_v_local[1])
            target_roll_deg = np.clip(target_roll_deg, -30, 30)

            # NOTE assume z to point downward, NED frame
            #print("target: %.2f, actual %.2f"%(target_v_local[2],self.drone_vel[this_id][2]))
            # NOTE using world frame z velocity
            target_thrust = self.baseThrust - self.vz_pids[this_id].control(target_v_local[2], self.drone_vel[this_id][2]) * self.thrustScale
            target_thrust = int(np.clip(target_thrust,self.minThrust,0xFFFF))
            target_yawrate_deg_s = 0
            cf.commander.send_setpoint(target_roll_deg,-target_pitch_deg,-target_yawrate_deg_s,target_thrust)
            self.log_dict['thrust'] = target_thrust
            #print(target_roll_deg,-target_pitch_deg,target_thrust)
            #print("thrust = %d"%target_thrust)
        except Exception as e:
            print_error(e)
        #print("\t roll: %.1f, pitch: %.1f, yawrate: %.1f, thrust %d"%(target_roll_deg,target_pitch_deg,target_yawrate_deg_s,target_thrust))
        return

    def wait_for_param_download(self,scf):
        while not scf.cf.param.is_updated:
            time.sleep(1.0)
        print('Parameters downloaded for', scf.cf.link_uri)

    # get current vehicle state
    # send to velocity planner
    # call velocity controller
    # update velocity commands
    def controlThread(self,):
        '''
        t0 = time()
        # parameters for a circle
        radius = [0.6,0.3]
        T = 10.0
        pi = np.pi
        '''

        while not self.quit_flag.isSet():
            ret = self.new_state.wait(0.1)
            if (not ret):
                continue
            self.new_state.clear()

            reduced_states = [[s[0],s[1]] for s in self.drone_states]
            # call planar controller (Ahmad's)
            velocity_commands_planar = self.multiAgentControl(reduced_states)

            # actuate 
            for i in range(self.drone_count):
                cmd = Planar(vx=velocity_commands_planar[i][0], vy=velocity_commands_planar[i][1], z=self.target_z_array[i])
                print(velocity_commands_planar[i][0],velocity_commands_planar[i][1])
                self.commands[i].put(cmd)
                self.new_command[i].set()
            '''
            # control a circle
            ts = time()
            for i in range(self.drone_count):
                target_x = radius[i] * cos( 2*pi/T*(ts - t0))
                target_y = radius[i] * sin( 2*pi/T*(ts - t0))
                cmd = Pos(x=target_x, y=target_y, z=self.target_z_array[0])
                self.commands[i].put(cmd)
                self.new_command[i].set()
            '''



if __name__ == '__main__':
    ins = MultipleDrone()
    ins.run()
