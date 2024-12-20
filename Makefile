GIT_HASH = $(shell git rev-parse --short HEAD)
GIT_TAG = $(shell git describe --tags --exact-match 2>/dev/null || echo "No Tag")

ifeq ($(GIT_TAG), No Tag)
	VERSION_ID = "commit-hash:$(GIT_HASH)"
else
	VERSION_ID = "$(GIT_TAG)"
endif

COMPILER = nvcc
COMPILER_FLAGS = -gencode arch=compute_89,code=sm_89 -DVERSION="\"$(VERSION_ID)\""
LINKER_FLAGS =

SOURCE_FILES = $(wildcard ./src/*.cu)
OBJECT_FILES = $(SOURCE_FILES:./src/%.cu=./obj/%.o)

TARGET_FILE = ./dust_collider

all: $(TARGET_FILE) clean
	@echo ""
	@echo "### COMPILATION COMPLETE"
	@echo "The executable is $(TARGET_FILE)"

$(TARGET_FILE): $(OBJECT_FILES)
	@echo $(OBJECT_FILES)
	@echo ""
	@echo "### LINKING OBJECT FILES"
	$(COMPILER) $(LINKER_FLAGS) $(OBJECT_FILES) -o $(TARGET_FILE)

./obj/%.o: ./src/%.cu
	@echo ""
	@echo "### COMPILING SOURCE FILE: $<"
	@mkdir -p ./obj
	$(COMPILER) $(COMPILER_FLAGS) -c $< -o $@

clean:
	@echo ""
	@echo "### REMOVING OBJECT DIRECTORY"
	rm -rf ./obj

very-clean: clean
	@echo ""
	@echo "### REMOVING EXECUTABLE"
	rm -f $(TARGET_FILE)
