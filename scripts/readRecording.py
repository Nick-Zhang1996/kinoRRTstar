import pickle
import matplotlib.pyplot as plt
import numpy as np

f = open("record.p","rb")
data = pickle.load(f)
f.close()

data = np.array(data)
x = data[:,0,0]
y = data[:,0,1]
z = data[:,0,2]
rx = data[:,0,3]
ry = data[:,0,4]
rz = data[:,0,5]

# Plot 3d.
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z, color='r', label='actual')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('altitude (z)')
ax.legend()
plt.show()
