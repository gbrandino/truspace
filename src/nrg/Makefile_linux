ROOT=${PWD}/../..
INSTALL_DIR=${ROOT}/lib
HEADER_DIR=${ROOT}/include
BLAS=${ROOT}/dep/OpenBLAS-0.3.17
BLAS_LIB=-lopenblas
LAPACK=${ROOT}/dep/OpenBLAS-0.3.17# it is better not to have this field empty, to prevent the -rpath to use /lib. If your blas distribution contains lapack, make it point to openblas root folder
LAPACK_LIB=
SCALAPACK_LIB=${ROOT}/scalapack-2.0.0/libscalapack.a
MPICXX= mpicxx

CFLAGS=-I${BLAS}/include -std=c++11 -fopenmp -O2 -fPIC
LDFLAGS = -shared 
LIBS= -L${BLAS}/lib/  ${BLAS_LIB}  -L${LAPACK}/lib/ ${LAPACK_LIB}  -lpthread -lstdc++ -fopenmp
SRC=parameters.cpp states.cpp matrix_elements.cpp filenames.cpp
SRC_SERIAL=NRG.cpp Hamiltonian.cpp
SRC_MPI=NRG_MPI.cpp Hamiltonian_MPI.cpp 

OBJS=$(SRC:.cpp=.o)
OBJS_SERIAL=$(SRC_SERIAL:.cpp=.o)
OBJS_MPI=$(SRC_MPI:.cpp=.o)

serial: libnrg.so

mpi : libnrg_mpi.so 

$(OBJS): %.o : %.cpp %.h config.h
	$(CXX) $(CFLAGS) -c $<

$(OBJS_SERIAL): %.o : %.cpp %.h config.h
	$(CXX) $(CFLAGS) -c $<


$(OBJS_MPI): %.o : %.cpp %.h config.h
	$(MPICXX) -DMPI_PAR $(CFLAGS) -c $<


libnrg.so: $(OBJS_SERIAL) $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^
	mkdir -p ${INSTALL_DIR}
	cp $@ ${INSTALL_DIR}
	mkdir -p ${HEADER_DIR}
	cp *.h ${HEADER_DIR} 

libnrg_mpi.so: $(OBJS_MPI) $(OBJS)
	$(MPICXX) $(LDFLAGS) -o $@ $^
	mkdir -p ${INSTALL_DIR}
	cp $@ ${INSTALL_DIR}
	mkdir -p ${HEADER_DIR}
	cp *.h ${HEADER_DIR}


clean:
	rm -f *.x *.o *.so *.dylib


