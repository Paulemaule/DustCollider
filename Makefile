GIT_BRANCH = branch
GIT_HASH = hash
GIT_TAG = tag

VERSION_ID = "$(GIT_BRANCH):$(GIT_HASH)"

COMPILER = g++
COMPILER_FLAGS = -O3 -DVERSION="\"$(VERSION_ID)\""
LINKER_FLAGS = 

SOURCE_FILES = $(wildcard ./src/*.cpp)
OBJECT_FILES = $(SOURCE_FILES:./src/%.cpp=./obj/%.o)

TARGET_FILE = ./dust_collider

all: $(TARGET_FILE) clean
	@echo ""
	@echo "### COMPILATION COMPLETE"
	@echo "The executable is $(TARGET_FILE)"

$(TARGET_FILE): $(OBJECT_FILES)
	@echo ""
	@echo "### LINKING OBJECT FILES"
	$(COMPILER) $(LINKER_FLAGS) $(OBJECT_FILES) -o $(TARGET_FILE)

./obj/%.o: ./src/%.cpp
	@echo ""
	@echo "### COMPILING SOURCE FILE: $<"
	@mkdir -p ./obj
	$(COMPILER) $(COMPILER_FLAGS) -c $< -o $@

clean:
	@echo "### REMOVING OBJECT DIRECTORY"
	rm -rf ./obj

very-clean: clean
	@echo "### REMOVING EXECUTABLE"
	rm -f $(TARGET_FILE)