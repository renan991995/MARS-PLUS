# MARS+: Molecular Assembly and Representation Suite - Plus



## The old MARS

In computer-aided molecular design (CAMD) [1][2], the capability of generating new molecular species from existing one is vital. MARS program [3] is devised for such demand. It consists of two components: 

* **Molecular data structure (MDS):** Base elements and 5 arrays of integers
* **Genetic operators:** ring formation, addition, subtraction, exchange, crossover, and combination

To initiate a MARS task, one should input the 3D structures of the starting molecules. These structures will then be converted into MDS representation, where a structure is recognized as a network of base elements. Each of the genetic operators will be applied to each of the possible substructures in each of the MDSs. As a result, a number of new species can be generated.  



## What's new in MARS+?

MARS+ is based on the prototype of MARS [3], with various improvements:

1. **The expansion of base element library**
        1-1. Group-like elements are allowed now.
        1-2. Common neutral atoms, ionic cores, and anionic cores are included.

2. **The generalization of MDS**
		2-1. An extra array of integers is used to bookkeep atomic chirality.
		2-2. Two extra arrays of integers are used to bookkeep cis-trans isomerism.
		2-3. An extra array of integers is used to bookkeep cyclic bonds.
		2-4. Multiple ring numbers on an atom are allowed now.
		2-5. The representation of 2-component chemical is allowed now. (1:1 ILs are demonstrated here)

3. **The generalization of genetic opertors**
		3-1. Refinement of old operators:
			 3-1-1. The feasibility of molecular connectivity is ensured after subtraction.
			 3-1-2. Multiple ring numbers on an atom can happen through cyclization operator.
			 3-1-3. Most of them are greatly revised using a more consistent approach.
			 3-1-4. The stability and flexibility are greatly enhanced. 
		3-2. Development of new operators: insertion, decyclization, element change, cis-trans inversion, chirality inversion, and component switch.
		3-3. With better reversibility for the new scheme of operators, undoing an operator is easier. This reduces the bias to certain types of molecules (e.g. polycyclics).
		3-4. Check cis-trans and chirality after genetic operations. (default: trans and clockwise winding)
		3-5. For imine substructure, indicate the lone pair of nitrogen atom by null atom "*".

    ![uni-molecular operations](./imgs/uni-molecular_operations.png "uni-molecular operations")
    
    ![crossover operation](./imgs/crossover.png "crossover")
    
    ![combination operation](./imgs/combination.png "combination")

4. **Wrapping Open Babel [4] functions into MARS+**
		4-1. This facilitates the inputting of starting structures to the program. Now one only needs to input SMILES.
		4-2. The perception for connectivity of inputted molecules are more robust.



## Development environment

* Linux CentOS 7
* GNU g++ compiler v9.2.0  (C++11)
* Open Babel v3.1.0 
* Eigen v3.3.7
* Cmake v3.15.5 
* Make v4.2 



## Compiling 

The MARS+ source code consists of 7 header files and 7 cpp files: (see `./src/` directory)

    ELEMENTS.h    MOLECULE.h    CASES_NEU.h    CASES_IL_INDEPENDENT.h    CASES_IL.h    UTILITY.h    PARAMETER.h
    ELEMENTS.cpp  MOLECULE.cpp  CASES_NEU.cpp  CASES_IL_INDEPENDENT.cpp  CASES_IL.cpp  UTILITY.cpp  main.cpp 
    
Before compiling, relevant softwares should be installed. For simplicity, it is recommended to get them through [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/download). If either of them is available, one can directly creates `MARS+` environment and imports the necessary packages from the provided YML file.

    conda env create --file ./MARS+_env.yml
    
Now activate the `MARS+` environment and compile. 

    conda activate MARS+
    cd ./src
    rm -r ./cmake_install.cmake ./CMakeFiles/ ./CMakeCache.txt 2> /dev/null
    cmake ./CMakeLists.txt
    make -j [N]                 ("-j [N]" is optional. Parallel compiling with N jobs at once)

An executable file `MARS-PLUS` will be generated in `./src/`.



## Usage

There are several input files for MARS+: (see `./INPUTS/` directory)
        
    ./INPUTS/control.in                     : controls the input, output, and calculation options.
    ./INPUTS/ELEMENT_LISTS/element_list.txt : a list that defines base element library.
    ./INPUTS/INPUT_CHEMICALS/IL4.txt        : the starting chemicals.
    
Please read the instructions in `./INPUTS/control.in` and `./INPUTS/ELEMENT_LISTS/element_list.txt`.
Make sure you have properly set the parameters. Now launch the MARS+:

	conda activate MARS+
	cd ./src/
	./MARS-PLUS ./INPUTS/control.in

Alternatively, you may use PBS scheduler if available.

	cd ./src/
	qsub ./job.sh

The results for each of the operations will be outputted as a file (see `./LOG_FILES/` directory).
For example, the results of bond change operation on an IL will be outputted to `./LOG_FILES/change_bnd_IL.txt`.



## Developers

Chen-Hsuan Huang (f07524028@ntu.edu.tw) and Shiang-Tai Lin (stlin@ntu.edu.tw).

Department of Chemical Engineering, National Taiwan University, Taipei, Taiwan



## References

[1] Austin, N. D.; Sahinidis, N. V.; Trahan, D. W., Computer-aided molecular design: An introduction and review of tools, applications, and solution techniques. Chem. Eng. Res. Des. 2016, 116, 2-26.

[2] Hsu, H. H.; Huang, C. H.; Lin, S. T., Fully Automated Molecular Design with Atomic Resolution for Desired Thermophysical Properties. Ind. Eng. Chem. Res. 2018, 57, (29), 9683-9692.

[3] Hsu, H.-H.; Huang, C.-H.; Lin, S.-T., New Data Structure for Computational Molecular Design with Atomic or Fragment Resolution. J. Chem. Inf. Model. 2019, 59, (9), 3703-3713.
(https://github.com/hsuhsuanhao/MARS)

[4] O’Boyle, N. M.; Banck, M.; James, C. A.; Morley, C.; Vandermeersch, T.; Hutchison, G. R., Open Babel: An open chemical toolbox. J. Cheminf. 2011. (https://github.com/openbabel/openbabel)





