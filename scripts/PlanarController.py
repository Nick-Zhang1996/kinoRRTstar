import numpy as np

# states: [[x0,y0],[x1,y1],[x2,y2]] for 3 drones
# output:
# target velocitvy for each drone [[vx0,vy0],[vx1,vy1],[vx2,vy2]] for 3 drones
# Example
def planarControllerExample(states):
    vx0 = 0
    vy0 = 0
    vx1 = 0
    vy1 = 0
    vx2 = 0
    vy2 = 0
    return [[vx0,vy0],[vx1,vy1],[vx2,vy2]]

def getBD():
    # define incidence matrix
    # define gain matris

    # for three agrents
    #D=np.matrix([[-1,0],[1,-1],[0,1]])
    #B = np.diag([0.1,0.1,0.1])*5

    # for two agents
    D=np.matrix([[-1],[1]])
    B = np.diag([1,0.5])*0.3
    return B,D

# TODO update this
def planarController(states):

    B,D = getBD()
    # converter state variable to numpy array
    X = np.asarray(states).reshape(4,1)

    # define orthogonal verctor conversion matrix
    S=np.matrix([[0,-1],[1,0]])

    # calculate Z = kron(D',I)X
    Z = np.kron(np.transpose(D),np.eye(2))
    Z = np.matmul(Z,X)

    # calculate control u = kron(BD,S)Z
    BD=np.matmul(B,D)
    BDS=np.kron(BD,S)
    u = np.matmul(BDS,Z)

    return np.reshape(u.transpose(),(2,2)).tolist()
