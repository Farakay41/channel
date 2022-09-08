# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 14:24:54 2022

@author: Test
"""

import numpy as np
import matplotlib.pyplot as plt


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
    
dudy = 0.5496
    
start = 10

end = 39
count = 0


list_time = np.arange(start, end+1, 1, dtype=int)
print (list_time)
for i in list_time:
   duldy = []
   uv = []
   count += 1

   duldy = read_file_double("Data/duldy.%i.out" % i)
   uv = read_file_double("Data/uvs.%i.out" % i)
   
   product_a = np.zeros([len(duldy)], dtype=np.double)
   
   product_b = np.zeros([len(duldy)], dtype=np.double)
   
   for j in range(len(duldy)):
      product_a[j] = duldy[j] * ( -1 )* uv[j] / 1000   
      
   for j in range(len(duldy)):
      product_b[j] = dudy * ( -1 ) * uv[j]
      
      
   
   with open('Data/duldy_uvs.%i.out' % i , mode='w') as dafile:
    product_a.astype('float64').tofile(dafile)
   
   with open('Data/dudy_uvs.%i.out' % i , mode='w') as dafile:
    product_b.astype('float64').tofile(dafile)
    
print(count)













