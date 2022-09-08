# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 13:16:19 2022

@author: Test
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import diff


def read_file_double(filename):

    dtype_a = np.dtype(np.double)

    # Read the whole binary file
    with open(filename, "rb") as binary_file:

        # read the date into an array of double
        data = np.fromfile(binary_file, dtype_a)

        return data


def read_file_int(filename):

    dtype_a = np.dtype(np.intc)

    # Read the whole binary file
    with open(filename, "rb") as binary_file:

        # read the date into an array of double
        data = np.fromfile(binary_file, dtype_a)

        return data
    
ul = np.zeros([], dtype=np.double)
us = np.zeros([], dtype=np.double)    
y = np.zeros([1071], dtype=np.double)



ul= read_file_double("ul_at_ys.39.out")
us = read_file_double("us.39.out")


ul = ul.reshape(1536, 1071)
us = us.reshape(1536, 1071)




plt.figure()
plt.plot(ul[156])
plt.plot(y)
plt.grid(linestyle='dotted')
plt.xlabel('timestep')
plt.ylabel('Large-scale fluctuation')
plt.xlim([0, 1071])

plt.savefig("large.png",dpi=800)



















