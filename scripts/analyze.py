import pickle
import matplotlib.pyplot as plt
import numpy as np
logFilename = "./log.p"
output = open(logFilename,'rb')
data = pickle.load(output)
output.close()

data = np.array(data)
print("log length %.1f seconds"%(data.shape[0]/119.88))

x = data[:,0]
y = data[:,1]
z = data[:,2]
rx = data[:,3]
ry = data[:,4]
rz = data[:,5]

'''
vx = data[:,6]
vy = data[:,7]
vz = data[:,8]

thrust = data[:,9]
target_vx = data[:,10]
target_vy = data[:,11]
target_vz = data[:,12]
'''

# Plot 3d.
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(-x, -y, -z, color='r', label='actual')
ax.set_xlabel('-x')
ax.set_ylabel('-y')
ax.set_zlabel('altitude (-z)')
ax.legend()
plt.show()
