CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -pthread

TARGET_RELEASE := thermo_sph
TARGET_DEBUG   := thermo_sph_debug
TARGET_TEST    := run_tests

SRC_DIR   := src
TEST_DIR  := tests
BUILD_DIR := build

INCLUDES  := -I $(SRC_DIR) \
             -I $(SRC_DIR)/objects \
             -I $(SRC_DIR)/phys \
             -I $(TEST_DIR)

ALL_SRCS  := $(shell find $(SRC_DIR) -type f -name "*.cpp")
LIB_SRCS  := $(filter-out $(SRC_DIR)/main.cpp, $(ALL_SRCS))
TEST_SRCS := $(shell find $(TEST_DIR) -type f -name "*.cpp" 2>/dev/null)

OBJS_RELEASE := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/release/%.o, $(ALL_SRCS))
OBJS_DEBUG   := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/debug/%.o, $(ALL_SRCS))
OBJS_LIB_DBG := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/debug/%.o, $(LIB_SRCS))
OBJS_TEST    := $(patsubst $(TEST_DIR)/%.cpp, $(BUILD_DIR)/test/%.o, $(TEST_SRCS))

all: release

release: CXXFLAGS += -O3 -flto -DNDEBUG
release: $(TARGET_RELEASE)

$(TARGET_RELEASE): $(OBJS_RELEASE)
	@echo "Linking Release Binary: $@..."
	$(CXX) $(CXXFLAGS) $(OBJS_RELEASE) -o $@
	@echo "Release build successful."

$(BUILD_DIR)/release/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# --- DEBUG BUILD ---
debug: CXXFLAGS += -g -O0 -DDEBUG -fsanitize=address,undefined
debug: $(TARGET_DEBUG)

$(TARGET_DEBUG): $(OBJS_DEBUG)
	@echo "Linking Debug Binary: $@..."
	$(CXX) $(CXXFLAGS) $(OBJS_DEBUG) -o $@
	@echo "Debug build successful."

$(BUILD_DIR)/debug/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# --- UNIT TESTING ---
test: debug $(TARGET_TEST)
	@echo "\nRunning Unit Tests..."
	./$(TARGET_TEST)

$(TARGET_TEST): CXXFLAGS += -g -O0 -DDEBUG -fsanitize=address,undefined
$(TARGET_TEST): $(OBJS_LIB_DBG) $(OBJS_TEST)
	@echo "Linking Test Runner: $@..."
	$(CXX) $(CXXFLAGS) $(OBJS_LIB_DBG) $(OBJS_TEST) -o $@

$(BUILD_DIR)/test/%.o: $(TEST_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# --- UTILITIES ---
run: release
	./$(TARGET_RELEASE)

clean:
	@echo "Cleaning build directory and binaries..."
	rm -rf $(BUILD_DIR) $(TARGET_RELEASE) $(TARGET_DEBUG) $(TARGET_TEST)

.PHONY: all release debug test run clean