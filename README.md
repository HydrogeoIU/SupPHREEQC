SupPHREEQC
-----------------------------------------------------


Introduction
--------------


SupPHREEQC is a program written in C++11, which links PHREEQC and SUPCRTBL programs to facilitate geochemical modeling at elevated temperatures and pressures suitable for clastic and carbonate diagenesis, geological carbon storage, and geothermal applications. Users could read the manuscript “SupPHREEQC: A program to generate customized PHREEQC thermodynamic databases from SUPCRTBL and extend calculations to elevated pressures and temperatures” (Zhang et al., in review) for more details. Below, the contents of the software package, compilation, and operation are discussed.


Contents of the Software Package 
---------------------------------------


Core C++11 source code files:

	AqueousModel.h	FileInput.h	FileOutput.h
	Global.h		ScreenInput.h	ScreenOutput.h
	SpeciesMineralGas.h	Utility.h 
	AqueousModel.cpp	FileInput.cpp	FileOutput.cpp
	Global.cpp		main.cpp	ScreenInput.cpp
	ScreenOutput.cpp	SpeciesMineralGas.cpp	Utility.cpp


C++ library files in folder ‘OtherCppLib’ for 2D scattered data linear interpolation , which are from John Burkardt (https://people.sc.fsu.edu/~jburkardt/cpp_src/pwl_interp_2d_scattered/pwl_interp_2d_scattered.html):

	pwl_interp_2d_scattered.hpp	r8lib.hpp
	pwl_interp_2d_scattered.cpp	r8lib.cpp


Necessary data files in folder ‘data’:

	abnormal gamma.txt
	addreplace.txt
	addreplace-PHREEQCDAT.txt
	bdot vs TP.txt
	rxn.dat
	SOLUTION_MASTER_SPECIES.txt
spronsbl.dat
vmcriticalpoint.txt
kineticsScript.txt


Compiling
-------------


GNU Scientific Library (GSL) is required when compiling. SupPHREEQC has been successfully compiled on Windows 10 and CentOS 7 using GCC.


Running
----------


cpronsbl.exe and SUPCRTBL.exe are required before running. Their source codes can be downloaded from https://github.com/HydrogeoIU/SUPCRTBL.


When running, please place SupPHREEQC.exe, cpronsbl.exe, SUPCRTBL.exe, and folder ‘data’ mentioned above in the same folder. Then run SupPHREEQC.exe and follow the software prompts. The expected PHREEQC thermodynamic databases can be created. For more detailed information, users can read the manuscript Zhang et al. (in review).


Questions Regarding the Software Package
----------------------------------------------------

Questions about the software package could be sent to:

	Guanru Zhang, Chengdu University of Technology, China, guanru.zhang@hotmail.com 
	Peng Lu, Saudi Aramco, Saudi Arabia, lvpeng00@gmail.com 
	Chen Zhu, Indiana University, USA, chenzhu@indiana.edu 


SupPHREEQC source codes and executable programs can be downloaded from https://scholarworks.iu.edu/dspace/handle/2022/23355. 
