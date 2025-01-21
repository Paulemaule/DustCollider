BUILD ?= Release

ifeq ($(BUILD), Release)
	BUILD_FLAGS = -DRELEASE -O3
else ifeq ($(BUILD), Debug)
	BUILD_FLAGS = -DDEBUG -g -G
else
	@echo This line is a hack and will cause make to crash. If this happens the BUILD version was unsupported.
endif

GIT_HASH = $(shell git rev-parse --short HEAD)
GIT_TAG = $(shell git describe --tags --exact-match 2>/dev/null || echo "No Tag")

ifeq ($(GIT_TAG), No Tag)
	VERSION_ID = "commit-hash:$(GIT_HASH)"
else
	VERSION_ID = "$(GIT_TAG)"
endif

SOURCE_DIR = ./src
BUILD_DIR = ./build

SOURCE_FILES = $(shell find $(SOURCE_DIR) -name '*.cu')
OBJECT_FILES = $(SOURCE_FILES:%=$(BUILD_DIR)/%.o)

DEPENDENCIES = $(OBJECT_FILES:.o=.d)
INCLUDE_DIRS = $(shell find $(SOURCE_DIR) -type d)
INCLUDE_FLAGS = $(addprefix -I,$(INCLUDE_DIRS)) -MMD -MP

TARGET_FILE = $(BUILD_DIR)/dust_collider

COMPILER = nvcc
COMPILER_FLAGS = -gencode arch=compute_89,code=sm_89 $(INCLUDE_FLAGS) -DVERSION="\"$(VERSION_ID)\"" $(BUILD_FLAGS)
LINKER_FLAGS =

all: $(TARGET_FILE)
	@echo "### COMPILATION COMPLETE"
	@echo "The executable is $(TARGET_FILE)"

$(TARGET_FILE): $(OBJECT_FILES)
	@echo "### LINKING OBJECT FILES"
	$(COMPILER) $(LINKER_FLAGS) $(OBJECT_FILES) -o $(TARGET_FILE)
	@echo ""

$(BUILD_DIR)/%.cu.o: %.cu
	@echo "### COMPILING SOURCE FILE: $<"
	@mkdir -p $(dir $@)
	$(COMPILER) $(COMPILER_FLAGS) -c $< -o $@
	@echo ""

-include $(DEPENDENCIES)

.PHONY: clean
clean:
	@echo "### REMOVING OBJECT DIRECTORY"
	rm -rf $(BUILD_DIR)
	@echo ""
