#!/bin/bash 

for Q in 7 ; do
    
    for LEVELS in 4; do

	# seq START STEP END
	for MASS in `seq -0.5 0.1 0.5`; do
	    ./example.sh ${MASS} ${LEVELS} ${Q}
	    ./deltalin ${MASS} ${LEVELS} ${Q}
	    
	done #MASS
    done #LEVELS
done #Q
