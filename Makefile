### Makefile - For compiling the dust collider code

## SETUP
# The Build configuration (Default: Release)
# BUILD can be either "Release" or "Debug"
BUILD ?= Release

# Determine additional compiler flags for different Builds
ifeq ($(BUILD), Release)
	BUILD_FLAGS = -DRELEASE
else ifeq ($(BUILD), Debug)
	BUILD_FLAGS = -DDEBUG -g -G
else
	@echo This line is a hack and will cause make to crash. If this happens the BUILD version was unsupported.
endif

# Git metadata retrieval
GIT_HASH = $(shell git rev-parse --short HEAD)
GIT_TAG = $(shell git describe --tags --exact-match 2>/dev/null || echo "No Tag")

# Pass git tag (or short commit hash, if no key is available) to compiler as VersionID
ifeq ($(GIT_TAG), No Tag)
	VERSION_ID = "commit-hash:$(GIT_HASH)"
else
	VERSION_ID = "$(GIT_TAG)"
endif

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
INCLUDE_FLAGS = $(addprefix -I,$(INCLUDE_DIRS)) -MMD -MP

# Path of output executable
TARGET_FILE = $(BUILD_DIR)/dust_collider

# Compiler configuration
COMPILER = nvcc
COMPILER_FLAGS = -gencode arch=compute_89,code=sm_89 $(INCLUDE_FLAGS) -DVERSION="\"$(VERSION_ID)\"" $(BUILD_FLAGS)
LINKER_FLAGS =

## TARGETS
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
	$(COMPILER) $(COMPILER_FLAGS) -c $< -o $@
	@echo ""

# Include dependencies for tracking
-include $(DEPENDENCIES)

# Removes all files produced by this Makefile
.PHONY: clean
clean:
	@echo "### REMOVING OBJECT DIRECTORY"
	rm -rf $(BUILD_DIR)
