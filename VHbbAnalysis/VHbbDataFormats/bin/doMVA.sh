#!/bin/sh

echo "Do you want to train the MVA?"
read mva

echo "Do you want to redo all the trees"
read trees

if [ "$trees" == "yes" ];
    then
    for i in MC_files/*.txt ; do
	NAME=`echo $i | perl -pe 's/.txt//'`
	echo $NAME
	make_tree $i $NAME >& $NAME.log &
	
    done
    
    for i in Data_files/*.txt ; do
	NAME=`echo $i | perl -pe 's/.txt//'`
	echo $NAME
	make_tree $i $NAME >& $NAME.log &
	
    done
fi

if [ "$mva"  == "yes" ];
    then

    JOBRUNNING=`ps -ef | grep make_BDT | grep -v grep | wc -l`
    echo $JOBRUNNING
    while [ $JOBRUNNING -gt 0 ]
      do
      echo "$JOBWORKING still running"
      echo "please wait..."
      sleep 60
      JOBRUNNING=`ps -ef | grep make_BDT | grep -v grep | wc -l`
    done

    echo "Which channel do you want to train? You can specify more than one. (Zmm, Zee, Znn, Wm, Wn, all)"
    read CHANNEL1 CHANNEL2 CHANNEL3 CHANNEL4 CHANNEL5

    if [ "$CHANNEL1" == "all" ]
	then
	echo "Running training for channel Zmm"
	TMVAClassification Zmm
	echo "Running training for channel Zee"
	TMVAClassification Zee
	echo "Running training for channel Znn"
	TMVAClassification Znn
	echo "Running training for channel Wm"
	TMVAClassification Wm
	echo "Running training for channel We"
	TMVAClassification We
    fi

    if [ "$CHANNEL1" ]
	then
	echo "Running training for channel $CHANNEL1"
	TMVAClassification $CHANNEL1
    fi
    if [ "$CHANNEL2" ]
	then
	echo "Running training for channel $CHANNEL2"
	TMVAClassification $CHANNEL2
    fi
    if [ "$CHANNEL3" ]
	then
	echo "Running training for channel $CHANNEL3"
	TMVAClassification $CHANNEL3
    fi
    if [ "$CHANNEL4" ]
	then
	echo "Running training for channel $CHANNEL4"
	TMVAClassification $CHANNEL4
    fi
    if [ "$CHANNEL5" ]
	then
	echo "Running training for channel $CHANNEL5"
	TMVAClassification $CHANNEL5
    fi
fi
