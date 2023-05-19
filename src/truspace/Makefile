ROOT=${PWD}/../..
CLN=${ROOT}/dep/cln-1.3.6
GINAC=${ROOT}/dep/ginac-1.8.6
BLAS=${ROOT}/dep/OpenBLAS-0.3.17
BLAS_LIB=-lopenblas
LAPACK=${ROOT}/dep/OpenBLAS-0.3.17# it is better not to have this field empty, to prevent the -rpath to use /lib. If your blas distribution contains lapack, make it point to openblas root folder
LAPACK_LIB=
SCALAPACK_LIB=${ROOT}/scalapack-2.0.0/libscalapack.a
NRG_LIB= -L${ROOT}/lib -lnrg
NRG_LIB_MPI=  -L${ROOT}/lib -lnrg_mpi
MPICXX= mpicxx

CFLAGS=-I${CLN}/include -I${GINAC}/include -I${BLAS}/include -I${ROOT}/include -std=c++11 -fopenmp -O2
LIBS= -L${GINAC}/lib/ -lginac -L${CLN}/lib/ -lcln -L${BLAS}/lib/  ${BLAS_LIB}  -L${LAPACK}/lib/ ${LAPACK_LIB} -lpthread -lstdc++ -fopenmp
SRC=characters.cpp gamma.cpp cpt.cpp block.cpp worker.cpp element_computation.cpp minimal_model_states.cpp diag.cpp minimal_model_filenames.cpp chirME.cpp minimal_model_parameters.cpp
SRC_FIELD= gamma.cpp cpt.cpp

OBJS=$(SRC:.cpp=.o)
#OBJS_MPI=$(SRC_MPI:.cpp=.o)
OBJS_FIELD=$(SRC_FIELD:.cpp=.o)

.PHONY: clean test test_mpi

serial: truspace.x  fields_and_structure_constants.x

mpi : truspace_mpi.x fields_and_structure_constants.x

$(OBJS): %.o : %.cpp %.h minimal_model_config.h
	$(CXX) $(CFLAGS) -c $<

minimal_model_Hamiltonian.o: minimal_model_Hamiltonian.cpp minimal_model_Hamiltonian.h minimal_model_config.h
	$(CXX) ${CFLAGS} -c $<

minimal_model_Hamiltonian_MPI.o: minimal_model_Hamiltonian_MPI.cpp minimal_model_Hamiltonian_MPI.h minimal_model_config.h
	$(CXX) ${CFLAGS} -c $<

truspace.x : truspace.cpp $(OBJS) minimal_model_Hamiltonian.o
	${CXX} ${CFLAGS}  $^ ${NRG_LIB} ${LIBS}  -o  $@  -Wl,-rpath ${CLN}/lib/ -Wl,-rpath ${GINAC}/lib/ -Wl,-rpath ${BLAS}/lib/ -Wl,-rpath ${LAPACK}/lib -Wl,-rpath ${ROOT}/lib/
	mkdir -p ${ROOT}/bin
	cp $@ ${ROOT}/bin/

truspace_mpi.x : CFLAGS += -DMPI_PAR 
truspace_mpi.x : CXX=${MPICXX}
truspace_mpi.x : truspace.cpp $(OBJS) minimal_model_Hamiltonian_MPI.o 
	${CXX} ${CFLAGS}  $^ ${NRG_LIB_MPI} ${LIBS}  ${SCALAPACK_LIB} -lgfortran -o  $@  -Wl,-rpath ${CLN}/lib/ -Wl,-rpath ${GINAC}/lib/ -Wl,-rpath ${BLAS}/lib/ -Wl,-rpath ${LAPACK}/lib -Wl,-rpath ${ROOT}/lib/ 
	mkdir -p ${ROOT}/bin
	cp $@ ${ROOT}/bin

fields_and_structure_constants.x : fields_and_structure_constants.cpp $(OBJS_FIELD) 
	$(CXX)  -I${ROOT}/include -std=c++11 -O2  $^ -o $@
	mkdir -p ${ROOT}/bin
	cp $@ ${ROOT}/bin/        

clean:
	rm -f *.x *.o

test:
	@rm ${ROOT}/test/Eigenvalues_p_3_q_4_minLev_0_maxLev_10_Per_2_c_0.09_Subspace_0_Mom_0 || true
	@rm ${ROOT}/test/FINAL_Eigenvalues_p_3_q_4_minLev_0_maxLev_10_Per_2_c_0.09_Subspace_0_Mom_0 || true
	@rm ${ROOT}/test/store* || true
	@cp truspace.x ${ROOT}/test/ && cd ${ROOT}/test/ && ./truspace.x input-ising-model input-ising-nrg > /dev/null
	@pytest  --junitxml results.xml  ${ROOT}/test/compare.py

test_mpi:
	@rm ${ROOT}/test/Eigenvalues_p_3_q_4_minLev_0_maxLev_10_Per_2_c_0.09_Subspace_0_Mom_0 || true
	@rm ${ROOT}/test/FINAL_Eigenvalues_p_3_q_4_minLev_0_maxLev_10_Per_2_c_0.09_Subspace_0_Mom_0 || true
	@rm ${ROOT}/test/store* || true
	@cp truspace_mpi.x ${ROOT}/test/ && cd ${ROOT}/test/ && mpirun -np 4 ./truspace_mpi.x input-ising-model input-ising-nrg > /dev/null
	@pytest  --junitxml results.xml  ${ROOT}/test/compare_mpi.py

