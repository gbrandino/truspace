List of files on which first clean up was done
===============================================

- characters.c -> characters.cpp and characters.h
- cpt.c -> cpt.cpp and cpt.h 
- block.c -> block.cpp and block.h
- diag.c -> diag.cpp and diag.h
- element_computation.c -> element_computation.cpp and element_computation.h
- gamma.cpp -> gamma.cpp and gamma.h
- states.c -> states.cpp and states.h (heavily modified. Object oriented design) 
- me.c -> chirME.h and chirME.cpp  (heavily modified. Object oriented design)
- diag.c -> diag.h and diag.cpp ( modified into a class, to handle mangling problems)

New files

- worker.cpp and worker.h   -> class that takes care of calling states calculation, chiral ME calculation, and call NRG 
- parameters.cpp and parameters.h -> class that takes care of input reading, calculating sizes of verma modules and structure constant, as well as analysing symmetry subspaces.
- tensors.cpp and tensors.h -> helper class to handle 2D,3D,4d,5D and 6D tensor, with automatic allocation and deallocation, bound check ensuring memory contiguity.
- filenames.h and filenames.cpp  -> helper class that hold output file names. 
- Hamiltonian.h and Hamiltonian.cpp --> Class that handles the computation of Hamiltonian matrix elements
- NRG.h and NRG.cpp --> model-agnostinc NRG and sweep manager class


to do list

- check Expecation values 

Documentation status

block.h                   done
characters.h              done
chirME.h                  done 
config.h                  done 
cpt.h                     done 
diag.h                    done 
element_computation.h     done   
filenames.h     
gamma.h                   done 
Hamiltonian.h       
NRG.h       
parameters.h              done 
states.h                  done
tensors.h                 done
worker.h                  done 
truspace.cpp              done

- create a section "expanding Truspace" with reference to doxygen files

NOTE:

- element computation seems to have some issues: valgrind reports a jump depending on an uninitialized value. However, data seems to be correct.
  If some problems arise, here is a good place to look
- element_computation.cpp requires "ulimit -s unlimited" to run for high truncation level ( not cool)

TO EXTEND SWEEP TO HIGHER MOMENTUM STATES

- modify states.cpp at line 840 ( the part of the code activated is nSweep > 0 )
- modify void Hamiltonian::sweep_block_BC ( line 192 of Hamiltonian.cpp) 
