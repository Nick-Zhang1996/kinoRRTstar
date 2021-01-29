# Tree class for RRT* management
from enum import Enum, auto
import numpy as np

# helper for indexing in Tree.nodes
class n(Enum):
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
    def __init__(self,root_state,max_node=2000):
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


    # add a node with prescribed node and children
    # return index of node
    def addNode(self,state,parent,children=[]):
        if (self.n_nodes == self.max_node):
            print_error("max node limit exceeded")
        self.nodes[self.n_nodes,:] = state
        self.parent[self.n_nodes] = parent
        self.children[self.n_nodes] = children
        self.n_nodes += 1
        return self.n_nodes-1

    # remove CHILD from NODE's children
    # NOTE the Tree class does not support removal of a node from tree
    # if a node is removed as a child here it must be rewired to another parent
    # we don't need old_parent since that can be deduced from child
    def transferChild(self,new_parent,child):
        old_parent = self.parent[child]
        self.nodes[old_parent].children.remove(child)
        # TODO should we maintain a sorted children list?
        self.nodes[new_parent].children.append(child)
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
    def getNeighbour(self,state,dist):
        xx = (self.nodes[:self.n_nodes,n.x] - state[0])**2
        yy = (self.nodes[:self.n_nodes,n.y] - state[1])**2
        zz = (self.nodes[:self.n_nodes,n.z] - state[2])**2
        dist_sq = (xx+yy+zz)
        mask = dist_sq < dist**2
        return mask.nonzero()



if __name__=="__main__":
    root = np.array((5,3,10,0,0,0,0),dtype=np.float)
    tree = Tree(root)
    
    state = np.array((4.5,3,10,0,0,0,0),dtype=np.float)
    n1 = tree.addNode(state,0)

    state = np.array((4,3.5,10,0,0,0,0),dtype=np.float)
    n2 = tree.addNode(state,n1)

    state = np.array((4.5,3.5,10,0,0,0,0),dtype=np.float)
    n3 = tree.addNode(state,n2)

    state = np.array((4.5,3.5,9,0,0,0,0),dtype=np.float)
    n4 = tree.addNode(state,n2)

    state = np.array((4,3,10,0,0,0),dtype=np.float)
    neighbour = tree.getNeighbour(state,0.6)
    print(neighbour)

    print(tree.getChildren(n2))
