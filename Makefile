# Compiler
CXX = g++

# Directories
SRC_DIR = src
BUILD_DIR = build

# Source files
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)

# Object files (replace src/ with build/ and .cpp with .o)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC_FILES))

# Flags
PYTHON ?= .venv/bin/python3
CXXFLAGS = -std=c++20 -Ofast -pthread -fPIC
CXXFLAGS += $(shell $(PYTHON) -m pybind11 --includes)

# Targets
all: $(OBJ_FILES)

# Build shared library (macOS)
lib: $(OBJ_FILES_NO_INTERFACE_NO_TEST)
	$(CXX) $(CXXFLAGS) -shared -o $(BUILD_DIR)/quboslib.so $(OBJ_FILES_NO_INTERFACE_NO_TEST)

# Rule to build object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

OBJ_FILES_NO_INTERFACE_NO_TEST := $(filter-out $(BUILD_DIR)/interface.o $(BUILD_DIR)/testing.o, $(OBJ_FILES))
# Core objects exclude interface.o, testing.o, and python_module.o (pybind11 bindings)
OBJ_FILES_CORE := $(filter-out $(BUILD_DIR)/interface.o $(BUILD_DIR)/testing.o $(BUILD_DIR)/python_module.o, $(OBJ_FILES))

# Compile interface.cpp into an executable, linking only core objects (no pybind11)
interface: $(BUILD_DIR)/interface.o $(OBJ_FILES_CORE)
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/interface $(BUILD_DIR)/interface.o $(OBJ_FILES_CORE)

# Compile testing.cpp into an executable, linking only core objects (no pybind11)
testing: $(BUILD_DIR)/testing.o $(OBJ_FILES_CORE)
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/testing $(BUILD_DIR)/testing.o $(OBJ_FILES_CORE)

# Clean build directory
clean:
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/interface $(BUILD_DIR)/testing $(BUILD_DIR)/libqubos.dylib
