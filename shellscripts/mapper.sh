#--------------------------------------------------------------------------------
#Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler
#
#Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
#released under the MIT licence.
#
#All right reserved.
#----------------------------------------------------------------------------------

#!/bin/bash

#Checking for existance for outout directory
if [ ! -d "../output" ]; then 
    mkdir ../output
fi

../build/qx_mapper ../layouts/ibmq_montreal.lut ../example/benchmark.qasm ../output/benchmark.qasm  1
