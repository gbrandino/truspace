ROOT=${PWD}/..
CLN=${ROOT}/dep/cln-1.3.4
GINAC=${ROOT}/dep/ginac-1.7.1
BLAS=${ROOT}/dep/OpenBLAS-0.3.17
BLAS_LIB=-lopenblas
LAPACK=${ROOT}/dep/OpenBLAS-0.3.17# it is better not to have this field empty, to prevent the -rpath to use /lib. If your blas distribution contains lapack, make it point to openblas root folder
LAPACK_LIB=
SCALAPACK_LIB=${ROOT}/scalapack-2.0.0/libscalapack.a
MPICXX= mpicxx

CFLAGS=-I${CLN}/include -I${GINAC}/include -I${BLAS}/include -std=c++11 -fopenmp -O2
CFLAGS_NRG= -fPIC
LDFLAGS = -shared 
LIBS= -L${GINAC}/lib/ -lginac -L${CLN}/lib/ -lcln -L${BLAS}/lib/  ${BLAS_LIB}  -L${LAPACK}/lib/ ${LAPACK_LIB}  -lpthread -lstdc++ -fopenmp
SRC=characters.cpp cpt.cpp gamma.cpp block.cpp worker.cpp element_computation.cpp states.cpp diag.cpp filenames.cpp chirME.cpp Hamiltonian.cpp minimal_model_parameters.cpp
SRC_MPI=characters.cpp cpt.cpp gamma.cpp block.cpp worker.cpp element_computation.cpp states.cpp diag.cpp filenames.cpp chirME.cpp Hamiltonian_MPI.cpp minimal_model_parameters.cpp
SRC_NRG=parameters.cpp NRG.cpp  
SRC_NRG_MPI=_parameters.cpp NRG_MPI.cpp 
SRC_FIELD= gamma.cpp cpt.cpp

OBJS=$(SRC:.cpp=.o)
OBJS_MPI=$(SRC_MPI:.cpp=.o)
OBJS_NRG=$(SRC_NRG:.cpp=.o)
OBJS_NRG_MPI=$(SRC_NRG_MPI:.cpp=.o)
OBJS_FIELD=$(SRC_FIELD:.cpp=.o)

all : libnrg.so truspace.x  fields_and_structure_constants.x

mpi : libnrg_mpi.so truspace_mpi.x fields_and_structure_constants.x

$(OBJS): %.o : %.cpp %.h config.h
	$(CXX) $(CFLAGS) -c $<

Hamiltonian_MPI.o : Hamiltonian_MPI.cpp Hamiltonian_MPI.h config.h
	$(CXX) $(CFLAGS) -c $<

NRG_MPI.o : NRG_MPI.cpp NRG_MPI.h config.h
	$(CXX) $(CFLAGS) $(CFLAGS_NRG) -c $<	

NRG.o : NRG.cpp NRG.h config.h
	$(CXX) $(CFLAGS) $(CFLAGS_NRG) -c $<

parameters.o : parameters.cpp parameters.h config.h
	$(CXX) $(CFLAGS) $(CFLAGS_NRG) -c $< 


libnrg.so: $(OBJS_NRG)
	$(CXX) $(LDFLAGS) -o $@ $^

libnrg_mpi.so: $(OBJS_NRG_MPI)
	$(CXX) $(LDFLAGS) -o $@ $^

truspace.x : truspace.cpp $(OBJS) libnrg.so
	${CXX} ${CFLAGS}  $^ ${LIBS} -o  $@  -Wl,-rpath ${CLN}/lib/ -Wl,-rpath ${GINAC}/lib/ -Wl,-rpath ${BLAS}/lib/ -Wl,-rpath ${LAPACK}/lib

truspace_mpi.x : CFLAGS += -DMPI_PAR 
truspace_mpi.x : CXX=${MPICXX}
truspace_mpi.x : truspace.cpp $(OBJS_MPI) libnrg_mpi.so
	${MPICXX} -DMPI_PAR ${CFLAGS}  $^ ${LIBS} ${SCALAPACK_LIB} -lgfortran -o  $@  -Wl,-rpath ${CLN}/lib/ -Wl,-rpath ${GINAC}/lib/ -Wl,-rpath ${BLAS}/lib/ -Wl,-rpath ${LAPACK}/lib 

fields_and_structure_constants.x : fields_and_structure_constants.cpp $(OBJS_FIELD) 
	$(CXX) -std=c++11 -O2  $^ -o $@
          

clean:
	rm -f *.x *.o


