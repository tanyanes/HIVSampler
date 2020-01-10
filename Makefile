#necessary C++ commands as compilation flags for warnings and debugging info
CXX=icpc
CXXFLAGS= -Wall -g

#paths to relevant directories that are being used by my program
INCLUDE = -I./include/ -I./gc -I./lib/cyCode
SOURCE = ./src
OBJ = ./obj
BIN = ./bin
DATA = ./data
LIB = ./lib

#paths to libraries and library flags
LDFLAGS = -mkl -qopenmp -pg -std=c++11 \
	#-lqhull

#define object files and source files
SRCS = $(SOURCE)/*.cpp
OBJS = $(SRCS:.cpp=.o)

#make sure that clean and test still run, even if files of the same name exist
.PHONY: clean test

#define variables for commands
EXEC = make

#define files being used

default: $(SRCS)
	 $(EXEC) -C src
	 $(EXEC) executable

executable: $(OBJ)/*.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) -o $(BIN)/executable $(OBJ)/*.o

clean:
	rm $(OBJ)/*.o
	rm $(BIN)/executable

test:
	#stuff for later
