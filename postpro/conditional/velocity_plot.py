# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 23:31:54 2022

@author: Test




"""
"""
import numpy as np
import matplotlib.pyplot as plt




def read_file(filename):

    dtype_a = np.dtype(np.double)

    # Read the whole binary file
    with open(filename, "rb") as binary_file:

        
        data = np.fromfile(binary_file, dtype_a)

        # reshape the data to fit the the four dimensions coordinates, positions x, y and z

    # I used information from the dns.in
        #data = data.reshape(3,96,129,67)
        data = data.reshape(8,60)
        return data


array_neg = np.empty([8
                      ,60], dtype= np.double)

array_neg = read_file("data.out")

print((array_neg[2][3]))

"""


import numpy as np
import matplotlib.pyplot as plt

#function reads binary file and writes to an array
def read_file(filename):

    dtype_a = np.dtype(np.double)

    # Read the whole binary file
    with open(filename, "rb") as binary_file:

        # read the date into an array of double
        data = np.fromfile(binary_file, dtype_a)

        return data



vel_pos = np.zeros([], dtype=np.double)
vel_neg = np.zeros([], dtype=np.double)

time_pos = np.zeros([], dtype=np.double)
time_neg = np.zeros([], dtype=np.double)


vel_pos = read_file("Velocity_pos_neg.out")
vel_neg= read_file("Velocity_neg_pos.out")


time_pos = read_file("time_pos_neg.out")
time_neg = read_file("time_neg_pos.out")


plt.figure()
plt.plot(time_pos, vel_pos)
plt.xlabel('tau')
plt.ylabel('Averaged_Velocity')
plt.title('Positive to Negative')
plt.show()


plt.figure()
plt.plot(time_neg, vel_neg)

plt.xlabel('tau')
plt.ylabel('Averaged_Velocity')
plt.title('Negative to Positive')

plt.show()






















