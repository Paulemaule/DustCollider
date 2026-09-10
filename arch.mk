### arch.mk - Shared CUDA architecture configuration
# This file contains code that ensures that the compiled code runs on any GPU 
# regardless of their compute capability.
# 
# Included by ./Makefile and ./tests/Makefile during compilation.
#
# The code is compiled into a single fat binary that holds machine code (SASS)
# for every architecture specified below. 
# The driver selects the matching one when the program starts, so one executable 
# runs on all of these GPUs without recompilation.
#
# Compile time and binary size scale with the length of the list.
# If the compute capability is known at compile time and fixed it can be specified via
#   make CUDA_ARCHS=89

CUDA_ARCHS ?= 80 86 89

# This ensures that the PTX of the newest architecture is also included.
PTX_ARCH := $(shell printf '%s\n' $(CUDA_ARCHS) | sort -n | tail -1)

# Build the compiler flags for:
# the Machine code for every requested architecture
ARCH_FLAGS := $(foreach arch,$(CUDA_ARCHS),-gencode arch=compute_$(arch),code=sm_$(arch))
# and the PTX of the newest one as a fallback for future GPUs.
ARCH_FLAGS += -gencode arch=compute_$(PTX_ARCH),code=compute_$(PTX_ARCH)
