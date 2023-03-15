#--------------------------------------------------------------------------------
#Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler
#
#Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
#released under the MIT licence.
#
#All right reserved.
#----------------------------------------------------------------------------------

import time, sys, os
import ntpath
import PyQMapping

PyQMapping.qx_mapper_file("../layouts/ibmq_montreal.lut", "../example/benchmark.qasm", "../output/benchmark2.qasm", 1) 
