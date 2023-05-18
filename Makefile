CC=/usr/local/cuda/bin/nvcc
ARCH=-c -O3 -arch sm_75 -std=c++11
MPIINC=-I/opt/openmpi/include -I/opt/openmpi/include/openmpi
INCLUDE=$(MPIINC) -I./include -I./include/algorithms
MPILIB=/usr/lib
BASE=src/algorithms/*.cu  src/kernels/*.cu src/util.cu src/cryptanalysis_types.cu src/automatic_linear_expansions.cu
FLAGS=-lmpi -lopen-rte -lopen-pal -lnsl -lutil -lm

BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj

SRC = $(wildcard $(BASE) src/crypt_dances_papers.cu src/crypt_dances_explorer.cu test/*.cu)
OBJ = $(SRC:%.cu=$(OBJ_DIR)/%.o)

.PHONY: all
all: $(BUILD_DIR)/crypt_dances_papers $(BUILD_DIR)/crypt_dances_explorer $(BUILD_DIR)/crypt_dances_tests

$(OBJ_DIR)/%.o: %.cu
	@mkdir -p $(@D)
	@$(CC) $(ARCH) $(INCLUDE) -rdc=true $< -o $@

$(BUILD_DIR)/crypt_dances_papers: $(OBJ)
	@mkdir -p $(BUILD_DIR)
	@$(CC) -lm -L$(MPILIB) $(FLAGS) $^ -o $@

$(BUILD_DIR)/crypt_dances_explorer: $(OBJ)
	@mkdir -p $(BUILD_DIR)
	@$(CC) -lm -L$(MPILIB) $(FLAGS) $^ -o $@

$(BUILD_DIR)/crypt_dances_tests: $(OBJ)
	@mkdir -p $(BUILD_DIR)
	@$(CC) -lm -L$(MPILIB) $(FLAGS) $^ -o $@

.PHONY : clean
clean :
	@rm -rf $(BUILD_DIR)
