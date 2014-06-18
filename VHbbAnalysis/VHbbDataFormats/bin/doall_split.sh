#!/bin/sh

for i in *.txt ; do
NAME=`echo $i | perl -pe 's/.txt//'`
echo $NAME
split -l 100 $i ${NAME}_split_ -d 

  for j in ${NAME}_split_* ; do
    make_histos $j $j >& $j.log &
  done

done
