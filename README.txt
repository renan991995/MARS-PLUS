########## Introduction to MARS ########## 
	MARS, Molecular Assembling and Representation Suite, [1] is a computer-aided molecular design (CAMD) [2]
	program for general purposes. The program uses five arrays of integers as the molecular data structure (MDS) 
	to bookkeep a molecular structure (i.e. constituent atoms, molecular connectivity, and formal charge etc.). 

	Genetic operators (i.e. ring formation, addition, subtraction, exchange, crossover, and combination) 
	were also developed so that a molecular data structure (MDS) can be used to modify in order to generate new 
	chemical speciess. MARS has been implemented in computer molecular design problems and has been found helpful. [3]

    Refs:
    [1] Hsu, H.-H.; Huang, C.-H.; Lin, S.-T., New Data Structure for Computational Molecular Design with Atomic or Fragment Resolution. 
		J. Chem. Inf. Model. 2019, 59, (9), 3703-3713.
        (https://github.com/hsuhsuanhao/MARS)

	[2] Austin, N. D.; Sahinidis, N. V.; Trahan, D. W., Computer-aided molecular design: An introduction and review of tools, applications, 
		and solution techniques. Chem. Eng. Res. Des. 2016, 116, 2-26.

    [3] Hsu, H. H.; Huang, C. H.; Lin, S. T., Fully Automated Molecular Design with Atomic Resolution for Desired Thermophysical Properties. 
		Ind. Eng. Chem. Res. 2018, 57, (29), 9683-9692.



########## Introduction to MARS-PLUS - What's new? ##########
	MARS-PLUS is a CAMD program for general purposes. 
    This program is developed based on the prototype of MARS [1], with various improvements:


	=========================================================================
	1. The expansion of base element library:
		1-1. Group-like elements are allowed now.
		1-2. Common neutral atoms, ionic cores, anionic cores are included.

	2. The generalization of MDS:
		2-1. An extra array of integers is used to bookkeep atomic chirality.
		2-2. Two extra arrays of integers are used to bookkeep cis-trans isomerism.
		2-3. An extra array of integers is used to bookkeep cyclic bonds.
		2-4. Multiple ring numbers on an atom are allowed now.
		2-5. The representation of 2-component chemical is allowed now. (1:1 ILs are demonstrated here)

	3. The improvements on genetic opertors:
		3-1. Refinement of old operators: 
			3-1-1. The feasibility of molecular connectivity is ensured after subtraction.
			3-1-2. Multiple ring numbers on an atom can happen through cyclization operator.
			3-1-3. Most of them are greatly revised for a more consistent approach.
		3-2. Development of new operators: 
			insertion, decyclization, element change, cis-trans inversion, 
			chirality inversion, and component switch.
		3-3. Check cis-trans and chirality after genetic operations. (default: trans and non-chiral)

    4. The incorporation of Open Babel facilitates the data utilization of designed chemicals.
		For example, one can convert a SMILES (outputted from MARS-PLUS) into 3D molecular structure 
		for ab initio calculations.
	
	5. Development of Open Babel-based SMILES enumerator:
		The SMILES enumerator can generate many synonymous SMILES from a given SMILES, and has been 
		found useful for improving NN-based models. [2] 
		We provided an Open Babel-based SMILES enumerator as a component of MARS-PLUS (see src/UTILITY.cpp).
		We also provided a standalone SMILES enumerator (see src/standalone_SMI_Enum/SMI_Enumerator.cpp) 
		that can compare the performance differences between Open Babel-based and RDKit-based enumerator.
	=========================================================================


    The source code consists of 7 header files and 7 cpp files: (see src/)
        ELEMENTS.h    MOLECULE.h    CASES_NEU.h    CASES_IL_INDEPENDENT.h    CASES_IL.h    UTILITY.h    PARAMETER.h
        ELEMENTS.cpp  MOLECULE.cpp  CASES_NEU.cpp  CASES_IL_INDEPENDENT.cpp  CASES_IL.cpp  UTILITY.cpp  main.cpp

    There are several input files to start MARS-PLUS: (see INPUTS/)
        INPUTS/control.in : controls the input, output, and calculation options.
        INPUTS/ELEMENT_LISTS/element_list.txt : a list that defines base element library.
        INPUTS/INPUT_CHEMICALS/IL4.txt : the beginning chemicals.

    The calculation results for each of the operations will be outputted as a file (see LOG_FILES/).
    For example, the results of bond change operation on an IL will be outputted to LOG_FILES/change_bnd_IL.txt


	Refs: 
    [1] Hsu, H.-H.; Huang, C.-H.; Lin, S.-T., New Data Structure for Computational Molecular Design with Atomic or Fragment Resolution.
        J. Chem. Inf. Model. 2019, 59, (9), 3703-3713.
        (https://github.com/hsuhsuanhao/MARS)

	[2] Bjerrum, E. J., SMILES Enumeration as Data Augmentation for Neural Network Modeling of Molecules. 2017.

			

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
	(*) RDKit 2020_03_1 (Q1 2020) Release
	(*) Boost v1.75.0		

	(*): It does not required for MARS-PLUS (src/UTILITY.cpp), 
		but is required for standalone SMILES enumerator (src/standalone_SMI_Enum/SMI_Enumerator.cpp).



########## Usage ##########
	1. MARS-PLUS
		Please read the instructions in INPUTS/control.in and INPUTS/ELEMENT_LISTS/element_list.txt .
		Make sure you have properly set the parameters before starting MARS-PLUS.

			cd src/
			rm -r ./cmake_install.cmake ./CMakeFiles/ ./CMakeCache.txt 2> /dev/null
			cmake ./CMakeLists.txt
			make
			./MARS-PLUS ./INPUTS/control.in



	2. Standalone SMILES enumerator
		Please read the instructions in src/standalone_SMI_Enum/control.in and src/standalone_SMI_Enum/SMI.txt .
		Make sure you have properly set the parameters before starting the standalone SMILES enumerator.

			cd src/standalone_SMI_Enum/
			rm -r ./cmake_install.cmake ./CMakeFiles/ ./CMakeCache.txt 2> /dev/null
			cmake ./CMakeLists.txt
			make
			./SMI_Enumerator ./control.in



	You might also utilize the job schedulers (e.g. PBS or Slurm) for both of the programs if available.



