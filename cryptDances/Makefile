#how to run:
#mpirun -np 17 -hostfile hosts ./papers_computational_results 

CC=/opt/cuda_11.1.1/bin/nvcc
ARCH=-c -O3 -arch sm_75
MPIINC=-I/opt/openmpi/include -I/opt/openmpi/include/openmpi
INCLUDE=$(MPIINC) -I./includes -I./includes/algorithms
MPILIB=/opt/openmpi/lib
BASE=src/algorithms/*.cu  src/kernels/*.cu src/util.cu src/cryptanalysis_types.cu src/automatic_linear_expansions.cu
FLAGS=-lmpi -lopen-rte -lopen-pal -lnsl -lutil -lm

.PHONY: all
all: crypt_dances_papers crypt_dances_explorer crypt_dances_tests

crypt_dances_papers: 
	$(CC) $(ARCH) $(INCLUDE) -rdc=true $(BASE) src/crypt_dances_papers.cu
	$(CC) -lm -L$(MPILIB) $(FLAGS) *.o -o crypt_dances_papers
	-rm $@.o
	
crypt_dances_explorer: 
	$(CC) $(ARCH) $(INCLUDE) -rdc=true $(BASE) src/crypt_dances_explorer.cu
	$(CC) -lm -L$(MPILIB) $(FLAGS) *.o -o crypt_dances_explorer
	-rm $@.o

crypt_dances_tests: 
	$(CC) $(ARCH) $(INCLUDE) -rdc=true $(BASE) tests/*.cu
	$(CC) -lm -L$(MPILIB) $(FLAGS) *.o -o crypt_dances_tests
	-rm $@.o

.PHONY : clean
clean :
	-rm -f crypt_dances_explorer crypt_dances_papers crypt_dances_tests *.o *.o
