import pickle
import matplotlib.pyplot as plt
import numpy as np
logFilename = "./log.p"
output = open(logFilename,'rb')
log_vec = pickle.load(output)
output.close()


x = []
y = []
z = []

rx = []
ry = []
rz = []

vx = []
vy = []
vz = []

thrust = []

log_vec = np.array(log_vec)
for i in range(log_vec.shape[0]):
    x.append( log_vec[i,:,0])
    y.append( log_vec[i,:,1])
    z.append( log_vec[i,:,2])

    rx.append( log_vec[i,:,3])
    ry.append( log_vec[i,:,4])
    rz.append( log_vec[i,:,5])

    vx.append( log_vec[i,:,6])
    vy.append( log_vec[i,:,7])
    vz.append( log_vec[i,:,8])

    thrust = log_vec[0,:,9]


#plt.plot(z)
#plt.plot(-thrust/200000)

for i in range(log_vec.shape[0]):
    #plt.plot(x[i],y[i])
    plt.plot(z[i])
plt.show()
