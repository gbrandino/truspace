ROOT=${PWD}/..
CLN=${ROOT}/dep/cln-1.3.4
GINAC=${ROOT}/dep/ginac-1.7.1
BLAS=${ROOT}/dep/OpenBLAS-0.2.19
BLAS_LIB=-lopenblas
LAPACK=${ROOT}/dep/OpenBLAS-0.2.19# it is better not to have this field empty, to prevent the -rpath to use /lib. If your blas distribution contains lapack, make it point to openblas root folder
LAPACK_LIB=
SCALAPACK_LIB=${ROOT}/scalapack-2.0.0/libscalapack.a
MPICXX= mpicxx

CFLAGS=-I${CLN}/include -I${GINAC}/include -I${BLAS}/include -std=c++11 -fopenmp -O2
LIBS= -L${GINAC}/lib/ -lginac -L${CLN}/lib/ -lcln -L${BLAS}/lib/  ${BLAS_LIB}  -L${LAPACK}/lib/ ${LAPACK_LIB}  -lpthread -lstdc++ -fopenmp
SRC= parameters.cpp characters.cpp cpt.cpp gamma.cpp block.cpp worker.cpp element_computation.cpp states.cpp diag.cpp filenames.cpp chirME.cpp Hamiltonian.cpp NRG.cpp
SRC_MPI=parameters.cpp characters.cpp cpt.cpp gamma.cpp block.cpp worker.cpp element_computation.cpp states.cpp diag.cpp filenames.cpp chirME.cpp Hamiltonian_MPI.cpp NRG_MPI.cpp
OBJS=$(SRC:.cpp=.o)
OBJS_MPI=$(SRC_MPI:.cpp=.o)

all : truspace.x

mpi : truspace_mpi.x

$(OBJS): %.o : %.cpp %.h config.h
	$(CXX) $(CFLAGS) -c $<

Hamiltonian_MPI.o : Hamiltonian_MPI.cpp Hamiltonian_MPI.h config.h
	$(CXX) $(CFLAGS) -c $<

NRG_MPI.o : NRG_MPI.cpp NRG_MPI.h config.h
	$(CXX) $(CFLAGS) -c $<	

truspace.x : truspace.cpp $(OBJS)
	${CXX} ${CFLAGS}  $^ ${LIBS} -o  $@  -Wl,-rpath ${CLN}/lib/ -Wl,-rpath ${GINAC}/lib/ -Wl,-rpath ${BLAS}/lib/ -Wl,-rpath ${LAPACK}/lib

truspace_mpi.x : CFLAGS += -DMPI_PAR 
truspace_mpi.x : CXX=${MPICXX}
truspace_mpi.x : truspace.cpp $(OBJS_MPI)
	${MPICXX} -DMPI_PAR ${CFLAGS}  $^ ${LIBS} ${SCALAPACK_LIB} -lgfortran -o  $@  -Wl,-rpath ${CLN}/lib/ -Wl,-rpath ${GINAC}/lib/ -Wl,-rpath ${BLAS}/lib/ -Wl,-rpath ${LAPACK}/lib 


clean:
	rm -f *.x *.o


