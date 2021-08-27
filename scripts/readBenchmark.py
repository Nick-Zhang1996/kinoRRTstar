import pickle
import numpy as np
import matplotlib.pyplot as plt

def show(tt, rrt_data, sst_data,  title, xlabel):
    plt.plot(tt, rrt_data, label='kinoRRT*')
    plt.plot(tt, sst_data, label='SST')
    plt.title(title)
    plt.xlabel("runtime (s)")
    plt.ylabel(xlabel)
    plt.legend()
    plt.show()

with open("log.p", 'rb') as f:
    data = pickle.load(f)
rrt, sst = data
rrt_nodes, rrt_cost, rrt_sols = rrt
sst_nodes, sst_cost, sst_sols = sst

rrt_cost[rrt_cost>1e3] = 0

tt = np.linspace(1,30,30)
show(tt, rrt_nodes, sst_nodes, "node count", "nodes")
show(tt, rrt_cost, sst_cost, "cost", "")
show(tt, rrt_sols, sst_sols, "solutions", "")
