(1) Introduction
		MARS-PLUS, Molecular Assembling and Representation Suite, is a program for
		general purpose computer aided molecular design. This program uses an 
		integer arry data structure to store the information of atom type, 
		atom connectivity, and bond order of a molecule. There are also 
		subroutines that can be used to modify the data structure in order to
		generate new molecules. This program is espicially useful for
		computer molecular design problems.

		The source code consists of 7 header files and 6 cpp files: (see src/)
			ELEMENTS.h		MOLECULE.h		CASES_NEU.h		CASES_IL_INDEPENDENT.h		CASES_IL.h		UTILITY.h		PARAMETER.h
			ELEMENTS.cpp	MOLECULE.cpp    CASES_NEU.cpp	CASES_IL_INDEPENDENT.cpp	CASES_IL.cpp	main.cpp

		There are several input files to start up MARS-PLUS: (see INPUTS/)
			INPUTS/control.in : 						controls the input, output, and calculation options.
			INPUTS/ELEMENT_LISTS/element_list.txt :		a list that defines base element library.
			INPUTS/INPUT_CHEMICALS/IL4.txt : 			the beginning chemicals.

		The calculation results for every operation will be outputted as a file (see LOG_FILES/):
			For example, the results of bond change operation on an IL will be outputted to LOG_FILES/change_bnd_IL.txt
			


(2) Development environment
		Linux CentOS 7
		g++ compiler from GNU Compiler Collection v9.2.0  	(or any compiler supporting C++11)
		Open Babel v3.1.0  									(compile from source code)
		cmake v3.15.5
		make v4.2



(3) Usage
		Please read the instructions in INPUTS/control.in and INPUTS/ELEMENT_LISTS/element_list.txt .
		Make sure you have properly set the parameters.

		cd src/
		rm -r ./cmake_install.cmake ./CMakeFiles/ ./CMakeCache.txt 2> /dev/null
		cmake ./CMakeLists.txt
		make
		./MARS ./INPUTS/control.in

		You might also utilize the job schedulers, such as Portable Batch System (PBS) and Simple Linux Utility for Resource Management (Slurm), if available.



(4) Developers
		This program is developed by Chen-Hsuan Huang and Shiang-Tai Lin.
		Computational Molecular Engineering Laboratory
		Department of Chemical Engineering, National Taiwan University, Taipei, Taiwan

		Correspondent: Shiang-Tai Lin (stlin@ntu.edu.tw)

