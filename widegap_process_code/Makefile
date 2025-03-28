
# Set the root directory
ROOT_DIR = ${ROOTSYS}

# Compiler settings
CXX = clang++
CC = clang

# Compiler flags
CXXFLAGS = -O3 -mtune=native `root-config --cflags` -Wall -fPIC -I. -stdlib=libc++

# Linker flags
LDFLAGS = -Wl,-rpath,${ROOT_DIR}/lib `root-config --libs` -O -stdlib=libc++

# Source files and object files
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
EXES = $(SRCS:.cpp=)

# Default target
all: $(EXES)

# Compile each source file and create the executable
%: %.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@

# Compile each source file into an object file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up object files and executables
clean:
	rm -f $(OBJS) $(EXES)

.PHONY: all clean
