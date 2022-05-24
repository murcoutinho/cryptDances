#how to run:
#mpirun -np 17 -hostfile hosts ./papers_computational_results 

CC=/opt/cuda_11.1.1/bin/nvcc
ARCH=-c -O3 -arch sm_75 -std=c++11
MPIINC=-I/opt/openmpi/include -I/opt/openmpi/include/openmpi
INCLUDE=$(MPIINC) -I./include -I./include/algorithms
MPILIB=/opt/openmpi/lib
BASE=src/algorithms/*.cu  src/kernels/*.cu src/util.cu src/cryptanalysis_types.cu src/automatic_linear_expansions.cu
FLAGS=-lmpi -lopen-rte -lopen-pal -lnsl -lutil -lm

BUILD_DIR = build

.PHONY: all
all: $(BUILD_DIR)/crypt_dances_papers $(BUILD_DIR)/crypt_dances_explorer $(BUILD_DIR)/crypt_dances_tests

$(BUILD_DIR)/crypt_dances_papers: 
	@mkdir -p $(BUILD_DIR)
	@$(CC) $(ARCH) $(INCLUDE) -rdc=true $(BASE) src/crypt_dances_papers.cu
	@$(CC) -lm -L$(MPILIB) $(FLAGS) *.o -o $@
	@rm *.o
	
$(BUILD_DIR)/crypt_dances_explorer: 
	@mkdir -p $(BUILD_DIR)
	@$(CC) $(ARCH) $(INCLUDE) -rdc=true $(BASE) src/crypt_dances_explorer.cu
	@$(CC) -lm -L$(MPILIB) $(FLAGS) *.o -o $@
	@rm *.o

$(BUILD_DIR)/crypt_dances_tests: 
	@mkdir -p $(BUILD_DIR)
	@$(CC) $(ARCH) $(INCLUDE) -rdc=true $(BASE) test/*.cu
	@$(CC) -lm -L$(MPILIB) $(FLAGS) *.o -o $@
	@rm *.o

.PHONY : clean
clean :
	@rm -rf $(BUILD_DIR)
	@rm -f *.o

