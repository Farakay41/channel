# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 12:19:14 2022

@author: Test
"""

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
        # data = data.reshape(3,96,129,67)
        data = data.reshape(8,60)
        return data


array_neg = np.empty([8
                      ,60], dtype= np.double)

array_neg = read_file("data.out")

print((array_neg[2][3]))

"""


# function reads binary file and writes to an array





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


n = 20
num = (2*n - 1)

start = 10

end = 39

vel_pos_ave = np.zeros([num], dtype=np.double)
vel_neg_ave = np.zeros([num], dtype=np.double)


time_pos = np.zeros([], dtype=np.double)
time_neg = np.zeros([], dtype=np.double)


list_time = np.arange(start, end+1, 1, dtype=int)
count = 0

print(list_time)

vel_pos = np.zeros([len(list_time), num], dtype=np.double)

vel_neg = np.zeros([len(list_time), num], dtype=np.double)

for i in list_time:

 vel_pos[i-start] = read_file_double("Velocity_pos_neg.%i.out" % i)
 vel_neg[i-start] = read_file_double("Velocity_neg_pos.%i.out" % i)

 count += 1

 for j in range(num):
      vel_pos_ave[j] += vel_pos[ (i-start), j ]
      vel_neg_ave[j] += vel_neg[(i-start), j ]
 
for j in range(num):

  vel_pos_ave[j] = vel_pos_ave[j] / len(list_time)
  print(vel_pos_ave[j])
  vel_neg_ave[j] = vel_neg_ave[j] / len(list_time)


array_pos = np.zeros([len(list_time)], dtype=np.double)
array_neg = np.zeros([len(list_time)], dtype=np.double)


std_array_pos = np.zeros([num], dtype=np.double)
std_array_neg = np.zeros([num], dtype=np.double)

for j in range(num):

  for i in list_time:
        array_pos[i-start] = vel_pos[(i-start), j]
        array_neg[i-start] = vel_neg[(i-start), j]

  std_array_pos[j] =  (np.std(array_pos)) / ( len(list_time) ** 0.5 )
 
  std_array_neg[j] =  (np.std(array_neg)) / ( len(list_time) ** 0.5 )
  

std_array_pos[:]  = 3 * std_array_pos[:] 
std_array_neg[:]  = 3 * std_array_neg[:] 


Result_pos = np.zeros([num, 2], dtype=np.double)
Result_neg  = np.zeros([num, 2], dtype=np.double)
  

Result_pos[:,0] = vel_pos_ave[:] - std_array_pos[:]
Result_pos[:,1] = vel_pos_ave[:] + std_array_pos[:]

Result_neg[:,0] = vel_neg_ave[:] - std_array_neg[:]
Result_neg[:,1] = vel_neg_ave[:] + std_array_neg[:]




print(Result_pos[:,0])
  
print(std_array_pos[3])  

print(std_array_neg[3]) 


time_pos = read_file_double("time_result.out")
time_neg = read_file_double("time_result.out")


t_pos= []

t_pos = time_pos[1:(num)]


t_neg= []

t_neg = time_neg[1:(num)]

print(len(t_pos))



dx = time_pos[5] - time_pos[4]

y = vel_pos_ave

j = vel_neg_ave
 

dy = np.zeros([], dtype=np.double)
dj = np.zeros([], dtype=np.double)


dy = diff(y)/dx
dj = diff(j)/dx



#np.savetxt('u_diff_P.txt', np.transpose([t_pos, dy]))

#np.savetxt('u_diff_N.txt', np.transpose([t_neg, dj]))



np.savetxt('duldy_p.txt', np.transpose([time_pos, vel_pos_ave, Result_pos[:,0], Result_pos[:, 1]]))

np.savetxt('duldy_n.txt', np.transpose([time_neg, vel_neg_ave, Result_neg[:,0], Result_neg[:, 1]]))



print(count)


print(len)

print(len (vel_pos))
print(len (time_pos))


print(len (time_neg))




# print(time_pos[7])
# print(vel_pos[37])

# print(vel_pos[67])

# print(vel_pos[97])


plt.figure()
plt.plot(time_pos, vel_pos_ave)

plt.xlabel('tau')
plt.ylabel('Averaged_Velocity_duldy_uvs')
# plt.xlim([-0.015, 0.015])
plt.title('Positive to Negative')

plt.figure()
plt.plot(time_neg, vel_neg_ave)

plt.xlabel('tau')
plt.ylabel('Averaged_Velocity_duldy_uvs')
# plt.xlim([-0.015, 0.015])
plt.title('Negative to Positive')


plt.figure()
plt.plot(time_pos, std_array_pos)

plt.xlabel('tau')
plt.ylabel('Standard deviation')
# plt.xlim([-0.015, 0.015])
#plt.title('Standard deviation')


plt.figure()
plt.plot(time_pos, std_array_neg)

plt.xlabel('tau')
plt.ylabel('Standard deviation')
# plt.xlim([-0.015, 0.015])
#plt.title('Standard deviation')



plt.show()






