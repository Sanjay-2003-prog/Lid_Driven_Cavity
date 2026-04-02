# Paths
AMREX_HOME := /media/sunder/Disk/sanjay/ImportCodes/amrex/build2d_cuda_hypre
HYPRE_DIR  := /media/sunder/Disk/sanjay/ImportCodes/hypre/src/hypre
#UMPIRE is required to run cuda with hypre
UMPIRE_DIR := /media/sunder/Disk/sanjay/ImportCodes/umpire_install


# Project directories
PROJ_DIR := .
SRC_DIR := $(PROJ_DIR)
BUILD_DIR := ./obj

# MPI paths
MPI_INC_DIR := /usr/lib/x86_64-linux-gnu/openmpi/include
MPI_LIB_DIR := /usr/lib/x86_64-linux-gnu/openmpi/lib

# Compilers
CXX = g++
NVCC = nvcc

# TINY_PROFILER flag 
TINY_PROFILE = TRUE

# CUDA architecture
CUDA_ARCH = sm_75

# C++ flags for regular compilation
CXXFLAGS = -std=c++17 -O3 -g -fPIC
CXXFLAGS += -DAMREX_USE_GPU -DAMREX_USE_CUDA -DAMREX_USE_MPI -DAMREX_USE_HYPRE
# Explicitly set architecture to your CPU's max supported features
CXXFLAGS += -march=cometlake -mtune=cometlake
#CXXFLAGS += -rdc=true #added for AMR

# ADD TINY_PROFILER FLAGS HERE
ifeq ($(TINY_PROFILE), TRUE)
    CXXFLAGS += -DAMREX_TINY_PROFILING
    NVCCFLAGS += -DAMREX_TINY_PROFILING
endif


# nvcc flags
NVCCFLAGS = -std=c++17 -O3 -g
NVCCFLAGS += -gencode arch=compute_75,code=sm_75 #changed due to SIGILL Error
#NVCCFLAGS += -arch=$(CUDA_ARCH)
NVCCFLAGS += --expt-relaxed-constexpr
NVCCFLAGS += --extended-lambda
NVCCFLAGS += -Xcompiler="-fPIC -march=x86-64 -mtune=generic"
NVCCFLAGS += -DAMREX_USE_GPU -DAMREX_USE_CUDA -DAMREX_USE_MPI
NVCCFLAGS += -rdc=true  # REQUIRED for device function linking

# Include paths
INCLUDES = -I. -I$(MPI_INC_DIR) -I$(AMREX_HOME)/include -I$(HYPRE_DIR)/include -I$(UMPIRE_DIR)/include

# Library paths
LIB_DIRS = $(AMREX_HOME)/lib $(MPI_LIB_DIR) /usr/local/cuda/lib64 $(HYPRE_DIR)/lib $(UMPIRE_DIR)/lib
#LIB_DIRS += $(AMREX_HOME)/lib -lamrex_gpu 
LIB_FLAGS = $(addprefix -L,$(LIB_DIRS))

# Libraries - ADD curand AND gfortran AND cuda and hypre
LIBS = -lamrex -lcurand -lcudart -lmpi -lgfortran -ldl -lrt -lpthread -lm 
LIBS += -lHYPRE -lcusparse -lcublas -lumpire -lcamp #these are required to run hypre

TARGET = run_cuda.exe

# Source files
SRCS = $(shell find $(SRC_DIR) -name "*.cpp")
OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)

# Default target
all: $(TARGET)

# Link target
$(TARGET): $(OBJS)
	$(NVCC) -o $@ $(OBJS) $(LIB_FLAGS) $(LIBS) $(NVCCFLAGS)

# Compile .cpp files: Use nvcc with -x cu flag
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(NVCC) -x cu -c -o $@ $(NVCCFLAGS) $(INCLUDES) $<

# Clean
clean:
	rm -f $(TARGET)
	rm -rf $(BUILD_DIR)

# Run
run: $(TARGET)
	mpiexec -n 1 ./$(TARGET)

.PHONY: all clean run
