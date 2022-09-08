#!/bin/bash

file_nos=$(seq 10 39)

for ii in ${file_nos}
do

  mv Data/ul_at_ys.${ii}.out ul.out
  mv Data/us.${ii}.out us.out
  #mv Data/vs.${ii}.out vs.out
  #mv Data/ws.${ii}.out ws.out
  #mv Data/uvs.${ii}.out uvs.out
  #mv Data/duldy.${ii}.out duldy.out
  #mv Data/duldy_uvs.${ii}.out duldy_uvs.out
  #mv Data/dudy_uvs.${ii}.out dudy_uvs.out

  mv Data/time.${ii}.out time.out
  g++ -o cut_1 cut_1.cpp -ljsoncpp
  ./cut_1
  mv Velocity_pos_neg.out Velocity_pos_neg.${ii}.out
  mv Velocity_neg_pos.out Velocity_neg_pos.${ii}.out
    


done