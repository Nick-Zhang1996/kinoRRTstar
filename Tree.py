# Tree class for RRT* management
from enum import IntEnum, auto
import numpy as np
from common import *
# DEBUG
#import copy

# helper for indexing in Tree.nodes
class n(IntEnum):
    x = 0
    y = 1
    z = 2
    vx = 3
    vy = 4
    vz = 5
    # cost start_node -> this_node
    cost = 6
    '''
    # parent index
    parent = 7
    # 1:true 0:false
    end_node = 8
    '''


# Tree class manages nodes
# each node has a unique index, assigned in ascending order when it's created, starting from 0
# the total number of nodes is Tree.n_nodes
# the index can be used in the following class vars
# self.nodes contains state and cost for each node
# self.parent contains a node's parent index parent_index = self.parent[child_index]
# self.children contains a node's children list, as a python list

# node with index 0 is root node

# NOTE
# self.nodes and other variables are preallocated, however, indexing beyond Tree.n_nodes causes undefined behavior
# because even though there may be state values pre-initialized for performance issues, the node have not be added to the tree (therefore invalid)
class Tree:
    def __init__(self,root_state,max_node=10):
        root_state = np.array(root_state, dtype=np.float).reshape(7)
        self.max_node = max_node
        self.nodes = np.zeros((max_node,7),dtype=np.float)
        # index of a parent
        self.parent = np.zeros(max_node,dtype= np.int)
        # 1 : true
        self.near_end_node = np.zeros(max_node,dtype= np.int)
        self.children = [[] for i in range(max_node)]
        self.n_nodes = 0

        self.addNode(root_state,0)
        return

    def getSolutionCount(self):
        return np.sum(near_end_node)

    # add a node with prescribed node and children
    # return index of node
    def addNode(self,state,parent):
        if (self.n_nodes == self.max_node):
            print_error("max node limit exceeded")
        self.nodes[self.n_nodes,:] = state
        self.parent[self.n_nodes] = parent
        #self.children[self.n_nodes] = []
        self.children[parent].append(self.n_nodes)
        self.n_nodes += 1
        return self.n_nodes-1

    # remove CHILD from NODE's children
    # NOTE the Tree class does not support removal of a node from tree
    # if a node is removed as a child here it must be rewired to another parent
    # we don't need old_parent since that can be deduced from child
    def transferChild(self,new_parent,child):
        old_parent = self.parent[child]
        self.children[old_parent].remove(child)
        # TODO should we maintain a sorted children list?
        self.children[new_parent].append(child)
        self.parent[child] = new_parent
        return

    # get list of indices of children
    # NOTE this does NOT create a copy, if modifications are done to the list
    # that may mess up the Tree class
    def getChildren(self,node):
        return self.children[node]

    # return number of total nodes
    def getNodeCount(self):
        return self.n_nodes

    # check if a node is connected to end node
    def isEnd(self,node):
        return self.near_end_node[node]

    # get all node index for a node within eucledian distance of dist
    # state can be any length 1d array, the first three values 
    # are expected to be x,y,z
    # return: neighbour indices
    def getNeighbour(self,state,dist):
        xx = (self.nodes[:self.n_nodes,n.x] - state[0])**2
        yy = (self.nodes[:self.n_nodes,n.y] - state[1])**2
        zz = (self.nodes[:self.n_nodes,n.z] - state[2])**2
        dist_sq = (xx+yy+zz)
        mask = dist_sq < dist**2
        return mask.nonzero()[0] 

    # find closest node index
    # and neighbours
    # similar to getNeighbour
    def findClosest(self,state,dist):
        xx = (self.nodes[:self.n_nodes,n.x] - state[0])**2
        yy = (self.nodes[:self.n_nodes,n.y] - state[1])**2
        zz = (self.nodes[:self.n_nodes,n.z] - state[2])**2
        dist_sq = (xx+yy+zz)
        index = np.argmin(dist_sq)

        mask = dist_sq < dist**2
        return index, mask.nonzero()[0] 

    # apply a delta in cost to node and all its descendents
    def spreadCostDelta(self,node,delta):
        # DEBUG
        #self.backup_nodes = self.nodes.copy()
        #self.backup_parent = self.parent.copy()
        #self.backup_children = copy.deepcopy(self.children)
        if (np.any(self.nodes[:,-1]<0)):
            print("error, negative cost")
            breakpoint()

        parent_queue = [node]

        count = 0
        while (len(parent_queue)>0):
            count += 1
            child_queue = [ child for parent in parent_queue for child in self.children[parent]]
            self.nodes[parent_queue,-1] += delta
            if (np.any(self.nodes[parent_queue,-1]) <0):
                breakpoint()

            parent_queue = child_queue
            if (count > 100):
                breakpoint()

        return

    def setNodeCost(self,nodecost):
        self.nodes[node,-1] = cost
        return

    # L2 distance from state (x,y,z,...) to node
    def dist2(self,state,node):
        x,y,z = self.nodes[node,:3]
        xx,yy,zz = state[:3]
        return ((x-xx)**2+(y-yy)**2+(z-zz)**2)**0.5

    # get UNIT vector from node to state (x,y,z,...)
    def vecFromNode(self,state,node):
        x,y,z = self.nodes[node,:3]
        xx,yy,zz = state[:3]
        vec = np.array([xx-x, yy-y, zz-z])
        dist = np.sum(vec**2)**0.5
        return vec/dist

    def getNodePos(self,node):
        return self.nodes[node,:3]

    def getNodeCost(self,node):
        return self.nodes[node,6]

    def getNodeState(self,node):
        return self.nodes[node,:6]

    def setEndState(self,node,state):
        self.near_end_node[node] = state
        return




if __name__=="__main__":
    root = np.array((5,3,10,0,0,0,0),dtype=np.float)
    tree = Tree(root)
    print(tree.children)
    
    state = np.array((4.5,3,10,0,0,0,0),dtype=np.float)
    n1 = tree.addNode(state,0)
    print(tree.children)

    state = np.array((4,3.5,10,0,0,0,0),dtype=np.float)
    n2 = tree.addNode(state,n1)
    print(tree.children)

    state = np.array((4.5,3.5,10,0,0,0,0),dtype=np.float)
    n3 = tree.addNode(state,n2)
    print(tree.children)

    state = np.array((4.5,3.5,9,0,0,0,0),dtype=np.float)
    n4 = tree.addNode(state,n2)
    print(tree.children)

    state = np.array((4,3,10,0,0,0),dtype=np.float)
    neighbour = tree.getNeighbour(state,0.6)
    print(neighbour)

    print(tree.getChildren(n2))
