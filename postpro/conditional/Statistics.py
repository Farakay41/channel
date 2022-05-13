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
array_pos = read_file("pos_neg.out")
array_neg = read_file("neg_pos.out")

#For the combined array
array = np.zeros([], dtype=np.intc)
array = np.append(array_pos, array_neg)

#sort from lowest to largest
array = np.sort(array, axis=None)



print('Kindly find the combined zero-crossing positions below: \n')

print(array)
print('\n')

#array to take in the intervals between successive zero_crossing positions
diff_array = np.zeros([len(array)-1], dtype=np.intc)

i = 1
while i < len(array):
    diff_array[i-1] = array[i] - array[i-1]
    i += 1
    
print('Kindly find the interval between succesive zero-crossing positions below: \n')

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

values = [value for value in range(min_diff, max_diff)]

probabilities = [dist.pdf(value) for value in values]

plt.hist(diff_array, bins=10, density=True)

plt.plot(values, probabilities)
plt.show()



model = KernelDensity(bandwidth=3, kernel='gaussian')

diff_array = diff_array.reshape((len(diff_array), 1))

model.fit(diff_array)


values = asarray([value for value in range (min_diff, max_diff) ])

values = values.reshape((len(values), 1))

probabilities = model.score_samples(values)
probabilities = exp(probabilities)

plt.hist(diff_array, bins=10, density=True)

plt.plot(values[:], probabilities)
plt.show()

