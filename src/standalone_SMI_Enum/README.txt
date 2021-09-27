########## Introduction to standalone SMILES enumerator ##########
	The SMILES enumerator can generate many synonymous SMILES from a given SMILES, and has been 
	found useful for improving NN-based models. [1] 
	We provided a standalone SMILES enumerator (see src/standalone_SMI_Enum/SMI_Enumerator.cpp) 
	that can compare the performance differences between Open Babel-based and RDKit-based enumerator.

    There are several input files to start the standalone SMILES enumerator: (see src/standalone_SMI_Enum/)
        src/standalone_SMI_Enum/control.in : controls the input, output, and calculation options.
        src/standalone_SMI_Enum/SMI.txt : the list of input SMILES'.

    The calculation results for each of the operations will be outputted as a file (see LOG_FILES/).
    For example, the results of bond change operation on an IL will be outputted to LOG_FILES/change_bnd_IL.txt


	Refs: 
	[1] Bjerrum, E. J., SMILES Enumeration as Data Augmentation for Neural Network Modeling of Molecules. 2017.

			

########## Developers ##########
    This program is developed by Chen-Hsuan Huang and Shiang-Tai Lin (stlin@ntu.edu.tw).
    Computational Molecular Engineering Laboratory
    Department of Chemical Engineering, National Taiwan University, Taipei, Taiwan



########## Development environment ##########
    Linux CentOS 7
    g++ compiler from GNU Compiler Collection v9.2.0        (or any compiler supporting C++11)
    Open Babel v3.1.0       (compile from source code)
    Cmake v3.15.5
    Make v4.2
	RDKit 2020_03_1 (Q1 2020) Release
	Boost v1.75.0		



########## Usage ##########
	1. Standalone SMILES enumerator
		Please read the instructions in src/standalone_SMI_Enum/control.in and src/standalone_SMI_Enum/SMI.txt .
		Make sure you have properly set the parameters before starting the standalone SMILES enumerator.

			cd src/standalone_SMI_Enum/
			rm -r ./cmake_install.cmake ./CMakeFiles/ ./CMakeCache.txt 2> /dev/null
			cmake ./CMakeLists.txt
			make
			./SMI_Enumerator ./control.in



	You might also utilize the job schedulers (e.g. PBS or Slurm) for both of the programs if available.



