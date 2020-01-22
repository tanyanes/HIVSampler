# necessary C++ commands as compilation flags for warnings and debugging info
CXX=icpc
CXXFLAGS= -Wall -g -O2 -ipo -xHost

# cyCode directories
CYCINC=-I../lib/cyCode

# QHULL Directories
QHULLDIR=/home/tanyanes/local
INCQHULLDIR= -I$(QHULLDIR)/include
LIBQHULLDIR=$(QHULLDIR)/lib
LDQHULL= -lqhull_r -lqhullcpp -L$(LIBQHULLDIR)
CQHULL=-I$(INCQHULLDIR)

#paths 
INCLUDE = -I../include/ $(CYCINC) $(INCQHULLDIR)
LDFLAGS= -mkl $(LDQHULL)
SOURCE= ../src
OBJ= ../obj
BIN= ../bin
DATA= ../data

SOURCES = ArithmeticHeader.cpp Stream.cpp newTest.cpp 
OBJS = newTest.o ArithmeticHeader.o Stream.o
EXEC = make

default: $(OBJS)
	mv *.o ../obj
	$(EXEC) executable

executable: $(OBJ)/*.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $(BIN)/executable $(OBJ)/*.o $(LDFLAGS)

newTest.o: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c newTest.cpp 

Stream.o: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c Stream.cpp -std=c++11

ArithmeticHeader.o: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c ArithmeticHeader.cpp 
