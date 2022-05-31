# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 14:23:31 2022

@author: Kayode
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import mean
from numpy import std
from scipy.stats import norm
import scipy.integrate as integrate
from numpy import hstack
from numpy import asarray
from numpy import exp

from sklearn.neighbors import KernelDensity
from sklearn.utils.fixes import parse_version

#function reads binary file and writes to an array
def read_file(filename):

    dtype_a = np.dtype(np.intc)

    # Read the whole binary file
    with open(filename, "rb") as binary_file:

        # read the date into an array of double
        data = np.fromfile(binary_file, dtype_a)

        return data



array_pos = np.zeros([], dtype=np.intc)
array_neg = np.zeros([], dtype=np.intc)

# writes the zero-crossing positions for both positive to negative and vice-versa

array_pos = read_file("ul.pnzc")
array_neg = read_file("ul.npzc")

print(len(array_pos) )



size_array_pos = read_file("ul_size.pnzc")
size_array_neg = read_file("ul_size.npzc")

print((size_array_pos[0:50]) )


#array= np.zeros([], dtype=np.intc)

diff_array = []
#diff_array = np.empty([], dtype=np.intc)

starting_neg = 0
starting_pos = 0



for i in range(len(size_array_pos)):
   
   if (size_array_pos[i] <= 1) :
     print(size_array_pos[i])
     starting_pos = starting_pos + size_array_pos[i]; 
     print(starting_pos)
 
        
   if (size_array_pos[i] > 1) :
       
     array_p = np.zeros([ size_array_pos[i] ], dtype=np.intc)
     #array_p =[]
     print(size_array_pos[i])
     print("sth")
     
     for j in range((size_array_pos[i])):
         
       array_p[j] =  array_pos[ starting_pos  + j ]
     #print(array_p)   
       
     starting_pos = starting_pos + size_array_pos[i]; 
     print(starting_pos)
   
   if  (size_array_neg[i] <= 1 ):
     print(size_array_neg[i])  
     starting_neg = starting_neg + size_array_neg[i]; 
     
    
   if (size_array_neg[i] > 1) :
     print(size_array_neg[i])  
     array_n = np.zeros([ size_array_neg[i] ], dtype=np.intc) 
     
     for j in range((size_array_neg[i])):
        array_n[j] =  array_neg[ starting_neg  + j ]
     #print(array_n)    
       # print(array_n)
     #print("f") 
     
     starting_neg = starting_neg + size_array_neg[i]; 
   
   #array= np.zeros([], dtype=np.intc)
   
   if (size_array_pos[i] + size_array_neg[i] >  2) :
       
     array = []
   
     if (size_array_pos[i] > 1) :             
      array = np.append(array, array_p)  
      print(array_p)
     if (size_array_neg[i] > 1) : 
      array = np.append(array, array_n) 
      print(array_n)
     array = np.sort(array, axis=None)
     print(array)
   
   
     diff_array_temp = np.zeros([(len(array)-1)], dtype=np.intc)
   
     for i in range((len(array) - 1)):
           
       diff_array_temp[i] = array[i+1] - array[i]
   #print(diff_array_temp) 
   
     diff_array = np.append(diff_array_temp, diff_array)

         



print('Kindly find the combined zero-crossing positions below: \n')

print(diff_array)
print('\n')




max_diff = np.max(diff_array)
print('The maximum interval is ' + str(max_diff) + '\n')

min_diff = np.min(diff_array)
print('The minimum interval is ' + str(min_diff) + '\n')

std_diff = np.std(diff_array)

average_diff = np.mean(diff_array) 
print('The average of the intervals is ' + str(round(average_diff, 3)) + '\n')

median_diff = np.median(diff_array)
print('The median of the intervals is ' + str(median_diff) + '\n')

First_quartile_diff = np.quantile(diff_array, 0.25)
print('The 1st quartile of the intervals is ' + str(First_quartile_diff) + '\n')

Third_quartile_diff = np.quantile(diff_array, 0.75)
print('The 3rd quartile of the intervals is ' + str(Third_quartile_diff) + '\n')

Int_Range_diff = Third_quartile_diff - First_quartile_diff

print('The interquartile range of the intervals is ' + str(Int_Range_diff) + '\n')



plt.hist(diff_array, bins=50)
plt.show()



dist = norm(average_diff, std_diff)

values = [value for value in range( int(min_diff), int(max_diff))]

probabilities = [dist.pdf(value) for value in values]

plt.hist(diff_array, bins=10, density=True)

plt.plot(values, probabilities)

plt.ylabel('PDF')
#plt.xlim([-0.015, 0.015])
plt.xlabel('Zero-crossings Intervals')
plt.title('PDF for Zero-crossings Intervals')

plt.show()



model = KernelDensity(bandwidth=3, kernel='gaussian')

diff_array = diff_array.reshape((len(diff_array), 1))

model.fit(diff_array)




values = asarray([value for value in range ( int(min_diff), int(max_diff)) ])

values = values.reshape((len(values), 1))

probabilities = model.score_samples(values)
probabilities = exp(probabilities)

plt.hist(diff_array, bins=10, density=True)

plt.plot(values[:], probabilities)
plt.show()

