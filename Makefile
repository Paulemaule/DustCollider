### Makefile - For compiling the dust collider code

## SETUP
# The Build configuration
# BUILD can be either "Release", "Debug" or "Test"
# If now BUILD is defined, default to Release
BUILD ?= Release

# Compiler flags per build type
BUILD_FLAGS.Release := -DRELEASE
BUILD_FLAGS.Debug   := -DDEBUG -g -G
BUILD_FLAGS.Test	:= -DTEST

BUILD_FLAGS := $(BUILD_FLAGS.$(BUILD))
ifeq ($(BUILD_FLAGS),)
	$(error Unsupported BUILD value '$(BUILD)'. Use Release or Debug.)
endif

# Helper script that defines the arch flags for GPU code compilation 
# to ensure compute capability compatibility.
# Defines 'ARCH_FLAGS'
include arch.mk

# Git version: tag if on a tag, tag+offset+hash if between tags, bare hash if no tags
VERSION_ID := $(shell git describe --tags --always 2>/dev/null)

# Directory structure
SOURCE_DIR = ./src
BUILD_DIR = ./build

# Determine Source and Object files
SOURCE_FILES = $(shell find $(SOURCE_DIR) -name '*.cu')
OBJECT_FILES = $(SOURCE_FILES:%=$(BUILD_DIR)/%.o)

# Dependency files for make's tracking
DEPENDENCIES = $(OBJECT_FILES:.o=.d)
# Source directories for proper header file search
INCLUDE_DIRS = $(shell find $(SOURCE_DIR) -type d)
# Additional compiler flags for header file search and dependency tracking
INCLUDE_FLAGS = $(addprefix -I,$(INCLUDE_DIRS))
DEP_FLAGS     = -MMD -MP

# Path of output executable
TARGET_FILE = $(BUILD_DIR)/dust_collider

# Compiler configuration
COMPILER = nvcc
COMPILER_FLAGS = -std=c++20 $(ARCH_FLAGS) $(INCLUDE_FLAGS) $(DEP_FLAGS) -DVERSION="\"$(VERSION_ID)\"" $(BUILD_FLAGS)
LINKER_FLAGS =

## TARGETS
.PHONY: all test clean

# Default target, builds the executable
all: $(TARGET_FILE)
	@echo "### COMPILATION COMPLETE"
	@echo "The executable is $(TARGET_FILE)"

# Linking the object files into a single executable
$(TARGET_FILE): $(OBJECT_FILES)
	@echo "### LINKING OBJECT FILES"
	$(COMPILER) $(LINKER_FLAGS) $(OBJECT_FILES) -o $(TARGET_FILE)
	@echo ""

# Compile all source files into object files
$(BUILD_DIR)/%.cu.o: %.cu
	@echo "### COMPILING SOURCE FILE: $<"
	@mkdir -p $(dir $@)
	$(COMPILER) $(COMPILER_FLAGS) -Xptxas -v -c $< -o $@
	@echo ""

# Include dependencies for tracking
-include $(DEPENDENCIES)

# Build and run the unit test suite
# Passes build and version information into test compilation.
test:
	$(MAKE) -C tests run BUILD=$(BUILD) VERSION_ID="$(VERSION_ID)"

# Removes all files produced by this Makefile (including test build artifacts)
clean:
	@echo "### REMOVING OBJECT DIRECTORY"
	rm -rf $(BUILD_DIR)
