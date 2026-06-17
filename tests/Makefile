### Makefile for the dust collider test suite

## SETUP
BUILD ?= Release

BUILD_FLAGS.Release := -DRELEASE
BUILD_FLAGS.Debug   := -DDEBUG -g -G
BUILD_FLAGS.Test	:= # No additional flag passed here, as the -DTEST if added as a fixed COMPILE_FLAG anyways

BUILD_FLAGS := $(BUILD_FLAGS.$(BUILD))
ifeq ($(BUILD_FLAGS),)
	$(error Unsupported BUILD value '$(BUILD)'. Use Release or Debug.)
endif

# Mirror the architecture flag from the main build.
ARCH_FLAGS = -gencode arch=compute_89,code=sm_89

# Include all subdirectories of src/ so test files can use the same
# #include paths as the main code (e.g. "utils/vector.cuh").
SRC_DIRS      := $(shell find ../src -type d)
INCLUDE_FLAGS = $(addprefix -I,$(SRC_DIRS)) -I.

COMPILER      = nvcc
COMPILE_FLAGS = $(ARCH_FLAGS) $(INCLUDE_FLAGS) $(BUILD_FLAGS) -DTEST

BUILD_DIR     = ../build/tests

## UNIT TEST TARGET
# All test_*.cu files are included into test_main.cu (single-TU build),
# so only test_main.cu needs to appear on the compiler command line.
UNIT_ENTRY  = test_main.cu
UNIT_DEPS   = $(wildcard test_definitions/unit/test_*.cu) framework/testkit.h
UNIT_TARGET = $(BUILD_DIR)/test_unit

## TARGETS
.PHONY: all unit run clean

all: unit

unit: $(UNIT_TARGET)
	@echo "### TEST BUILD COMPLETE: $(UNIT_TARGET)"

run: unit
	@echo "### RUNNING UNIT TESTS"
	@$(UNIT_TARGET); EXIT=$$?; \
	if [ $$EXIT -eq 0 ]; then echo "### ALL TESTS PASSED"; \
	else echo "### TESTS FAILED (exit $$EXIT)"; fi; \
	exit $$EXIT

$(UNIT_TARGET): $(UNIT_ENTRY) $(UNIT_DEPS)
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILE_FLAGS) $(UNIT_ENTRY) -o $(UNIT_TARGET)

clean:
	rm -rf $(BUILD_DIR)
