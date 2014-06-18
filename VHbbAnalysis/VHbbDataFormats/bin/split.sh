#!/bin/sh

for i in *.txt ; do
NAME=`echo $i | perl -pe 's/.txt//'`
echo $NAME
make_histos $i $NAME >& $NAME.log &

done
