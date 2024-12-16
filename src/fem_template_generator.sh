#!/bin/bash

maxcontact=$1

for (( i=0 ; i<$1 ; i++ ))
do
  sed s/THISID/$i/g $2/fem_impl_2.cpp > fem_impl_2_$i.cpp
done 





