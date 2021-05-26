# generate random number that we can reuse
import pickle
import numpy as np

f = open("random.p","wb")
data = []
for i in range(10000):
    new_pos = np.random.rand(3) * np.array([20,10,10])
    data.append(new_pos)
data = np.array(data)
print(data.shape)

pickle.dump(data,f)

