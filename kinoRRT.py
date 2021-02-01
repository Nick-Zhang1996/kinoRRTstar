# kinodynamic RRT star
import numpy as np
from math import log

from World import World
from DoubleIntegrator import DoubleIntegrator
from Tree import Tree,n
from timeUtil import execution_timer
import matplotlib.pyplot as plt

import pickle

debug_flag = False
def debug(*arg):
    global debug_flag
    if (debug_flag):
        print(*arg)
        #breakpoint()
    return

class kinoRRT:
    #
    # nodes: nodes to sample
    # if set to 0, quit after first solution is found
    def __init__(self, world, start_pos, end_pos, nodes=0, use_saved_rand = True):
        self.t = execution_timer(True)
        self.segment_len = 4.0

        #pos, vel, cost
        start_state = np.hstack([start_pos,[0,0,0],0]).astype(np.float)
        self.end_pos = np.array(end_pos,dtype=np.float)

        # init tree
        if (nodes == 0):
            self.tree = Tree(start_state,max_node=6000)
            self.max_nodes = 6000
            self.quit_after_first_solution = True
        else:
            self.tree = Tree(start_state,max_node=nodes)
            self.max_nodes = nodes
            self.quit_after_first_solution = False

        self.world = world
        self.di = DoubleIntegrator()

        # for determining best neighbour distance
        self.gamma = 40
        # problem dimension, 3D world: 3
        self.dim = 3
        self.min_neighbour_radius = 3.0

        # solutions contains information about nodes that can reach end_pos
        # each solution is an entry in self.solutions in for [node_id, total_cost_to_endnode]
        self.solutions = []
        self.found_solution = False

        self.use_saved_rand = use_saved_rand
        if (use_saved_rand):
            f = open("random.p","rb")
            self.saved_rand = pickle.load(f)
            self.saved_rand_index = 0

        return

    # run kinodynamic RRT for specified start and finish position
    # velocity at initial_pos is 0
    # velocity at final_pos is free
    # cost function is encoded in DoubleIntegrator class in the polynomials
    # and cannot be changed
    def run(self):
        for i in range(self.max_nodes-1):
            # randomly sample until we add a valid node to tree
            while (not self.sampleSpace()):
                #breakpoint()
                #print("bad node, resample")
                #print('.')
                pass
            #print("--- total nodes %d"%(self.tree.n_nodes))
            if (i % 100==0):
                print(i)


            # have we found a solution?
            if (self.found_solution and self.quit_after_first_solution):
                break

    # from original implementation
    def calcNeighbourDist(self):
        nun = self.tree.n_nodes
        ner = self.gamma * ( log(nun + 1.0) / nun )**(1.0 / self.dim);
        return min(ner,self.min_neighbour_radius)


    # sample space, if result is successful return True
    # if fail return False
    def sampleSpace(self):
        t = self.t
        t.s()

        t.s("sample")
        # randomly sample in world
        if (self.use_saved_rand):
            new_pos = self.saved_rand[self.saved_rand_index]
            self.saved_rand_index += 1
        else:
            new_pos = np.random.rand(3) * np.array(self.world.size)
        t.e("sample")

        debug("sampled point", new_pos)

        # find nearest point
        t.s("move 2 neighbour")
        closest_id, neighbour_id_vec = self.tree.findClosest(new_pos, self.calcNeighbourDist())
        # if too far move to that neighbour
        dist = self.tree.dist2(new_pos,closest_id)
        if ( dist > self.segment_len):
            debug("sampled point too far(%.3f), moving toward node %d"%(dist,closest_id))
            debug("node %d "%(closest_id)+str(self.tree.getNodePos(closest_id)))
            new_pos = self.tree.getNodePos(closest_id) + self.segment_len * self.tree.vecFromNode(new_pos,closest_id)
            debug("new pos "+str(new_pos)+"dist %.3f"%(self.tree.dist2(new_pos,closest_id)))

            neighbour_id_vec = self.tree.getNeighbour(new_pos, self.calcNeighbourDist())
            dist = self.segment_len
        t.e("move 2 neighbour")

        # collision?
        t.s("collision")
        if (not self.world.checkNoCollision(new_pos.reshape(1,3))):
            debug("collision with obstacle, resample...")
            # better luck next time
            t.e("collision")
            t.e()
            return False
        t.e("collision")

        # now we have a collision free point
        # and we know its nearest neighbour
        debug("no collision")

        # can we go from neighbour to new node?
        t.s("path collision")
        init_state = self.tree.getNodeState(closest_id)
        tf = self.di.DI3d_timeFreeVel(*init_state, *new_pos)
        state_vec = self.di.DI3d_stateFreeVel(tf,*init_state, *new_pos)
        # no need to check start and finish
        interior_state_vec = state_vec[1:-1]
        debug("checking path collision..")

        if (not self.world.checkNoCollision(interior_state_vec)):
            debug("path collision found")
            debug(interior_state_vec[:,:3])
            t.e("path collision")
            t.e()
            return False
        t.e("path collision")

        t.s("find parent")
        # no collision, add this sampled point to our tree
        # search neighbour for best parent
        debug("no path collision from %d"%(closest_id))
        best_parent = closest_id
        min_cost = self.di.DI3d_costFreeVel(tf, *init_state, *new_pos) + self.tree.getNodeCost(closest_id)
        assert(min_cost>0)
        debug("cost for new_node %.3f"%(min_cost))
        debug("checking neighbour for lower cost")
        debug("neighbours: "+str(neighbour_id_vec))

        # pos: pos
        # state: pos,vel
        # full_state: pos,vel,cost
        new_full_state = np.hstack([state_vec[-1],min_cost])
        for neighbour in neighbour_id_vec:
            debug("checking neighbour %d"%(neighbour))
            init_state = self.tree.getNodeState(neighbour)

            tf = self.di.DI3d_timeFreeVel(*init_state, *new_pos)
            state_vec = self.di.DI3d_stateFreeVel(tf,*init_state, *new_pos)
            # no need to check start and finish
            interior_state_vec = state_vec[1:-1]
            if (not self.world.checkNoCollision(interior_state_vec)):
                continue
            cost = self.di.DI3d_costFreeVel(tf, *init_state, *new_pos) + self.tree.getNodeCost(neighbour)
            assert(cost>0)
            debug("cost %.3f"%(cost))
            if (cost < min_cost):
                debug("better cost for %d %.3f < current min %.3f"%(neighbour, cost,min_cost))
                min_cost = cost
                best_parent = neighbour
                new_full_state = np.hstack([state_vec[-1],min_cost])
        t.e("find parent")

        # now add new_node to tree
        t.s("add")
        new_node_id = self.tree.addNode(new_full_state, best_parent)
        assert new_full_state[-1] > self.tree.getNodeCost(best_parent)
        new_state = new_full_state[:-1]
        debug("add new node(%d), parent %d, state"%(new_node_id,best_parent)+str(new_full_state))
        t.e("add")

        # is the new node close to end_state?
        t.s("end")
        end_pos = self.end_pos
        if (self.tree.dist2(end_pos,new_node_id) < self.segment_len):
            # can we reach end_state?
            tf = self.di.DI3d_timeFreeVel(*new_state, *end_pos)
            state_vec = self.di.DI3d_stateFreeVel(tf,*new_state, *end_pos)
            interior_state_vec = state_vec[1:-1]
            if (self.world.checkNoCollision(interior_state_vec)):
                # we've found a path
                self.tree.setEndState(new_node_id,1)
                cost = self.di.DI3d_costFreeVel(tf, *new_state, *self.end_pos) + self.tree.getNodeCost(new_node_id)
                assert(cost>0)

                self.solutions.append([new_node_id,cost])
                self.found_solution = True
                debug("new solution found with %d"%(new_node_id))
                # no need to rewire, new_node should remain a leaf node
                t.e("end")
                t.e()
                return True
        t.e("end")

        # rewire
        # can new node be a better parent for its neighbours? 
        t.s("rewire")
        debug("rewiring...")
        for neighbour in neighbour_id_vec:
            debug("checking possible child for new node :%d"%(neighbour))
            neighbour_state = self.tree.getNodeState(neighbour)
            tf = self.di.DI3d_time(*new_state, *neighbour_state)
            state_vec = self.di.DI3d_state(tf,*new_state, *neighbour_state)
            interior_state_vec = state_vec[1:-1]
            if (not self.world.checkNoCollision(interior_state_vec)):
                debug("collision")
                continue
            cost = self.di.DI3d_cost(tf, *new_state, *neighbour_state) + self.tree.getNodeCost(new_node_id)
            assert(cost>0)

            debug("no collision, cost %.3f"%(cost))
            if (cost < self.tree.getNodeCost(neighbour)):
                debug("new low cost %.3f"%(cost))
                d_cost = cost - self.tree.getNodeCost(neighbour)

                #print("transferChild(%d,%d)"%(new_node_id,neighbour))

                self.tree.transferChild(new_node_id,neighbour)
                self.tree.spreadCostDelta(neighbour,d_cost)
        t.e("rewire")
        debug("sample result in added node")
        t.e()
        return True


    # extract solutions from built tree
    # plot
    #   environment with obstacles
    #   sampled states
    #   trajectory through sampled states
    # if multiple solutions are found show each sequentially
    def showResult(self):
        sol = np.array(self.solutions)
        min_cost_id = np.argmin(sol[:,1])
        print("min cost %.2f"%(sol[min_cost_id,1]))

        # retrace steps to find the optimal trajectory
        final_index = self.solutions[min_cost_id][0]
        indices = [final_index]
        while (True):
            index = self.tree.parent[indices[-1]]
            indices.append(index)
            if (index == 0):
                break
        print("indices (reverse)"+str(indices))
        state_hist = []
        for i in range(len(indices)-1):
            init = self.tree.getNodeState(indices[len(indices)-1-i])
            final = self.tree.getNodeState(indices[len(indices)-2-i])

            tf = self.di.DI3d_time(*init, *final)
            state_vec = self.di.DI3d_state(tf,*init, *final)
            state_hist.append(state_vec[:-1,:])

        state_hist = np.array(state_hist).reshape(-1,6)
        pos_hist = state_hist[:,:3]

        node_pos = [self.tree.getNodePos(i) for i in indices[-1::-1]]
        node_pos = np.array(node_pos)

        ax = self.world.visualizeWorld(show=False)
        ax.plot(pos_hist[:,0], pos_hist[:,1], pos_hist[:,2],'b*')
        ax.scatter(node_pos[:,0], node_pos[:,1], node_pos[:,2], 'ro')
        plt.show()

        return


if __name__=="__main__":
    world = World(size=(20,10,10))
    obs = [[6,8],[0,5],[0,10]]
    world.addObstacle(obs)
    obs = [[6,8],[5,10],[0,5]]
    world.addObstacle(obs)
    obs = [[12,14],[5,10],[0,10]]
    world.addObstacle(obs)
    obs = [[12,14],[0,5],[5,10]]
    world.addObstacle(obs)
    #world.visualizeWorld()

    start_pos = (2,2,2)
    end_pos = (18,8,8)
    # keep searching until finding a solution
    rrt = kinoRRT(world,start_pos, end_pos,2000)
    rrt.run()
    rrt.showResult()
    print(rrt.solutions)
    rrt.t.summary()
