CXX=icpc
CXXFLAGS= -Wall -g


CYCINC=-I../lib/cyCode

# QHULL Directories
QHULLDIR=../qhull
INCQHULLDIR=$(QHULLDIR)/include
LIBQHULLDIR=$(QHULLDIR)/lib
LDQHULL=-lqhull -L$(LIBQHULLDIR)
CQHULL=-I$(INCQHULLDIR)


# CyCode Directories


LDQHUL=-lqhull -L$(QHULLLIBDIR)
#paths 
INCLUDE= -I../include -I../gc $(CYCINC)
LDFLAGS= -L./ -mkl
SOURCE= ../src
OBJ= ../obj
BIN= ../bin
DATA= ../data

SOURCES = ArithmeticHeader.cpp Stream.cpp newTest.cpp
OBJS = newTest.o ArithmeticHeader.o Stream.o

default: $(OBJS)
	mv *.o ../obj

newTest.o: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) -c newTest.cpp

Stream.o: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) -c Stream.cpp -std=c++11

ArithmeticHeader.o: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) -c ArithmeticHeader.cpp

#LFQHULL=-L/home/tanyanes/local/lib -lqhull
#LDFLAGS=$(LFQHULL)
#CFLAGS=-I/home/tanyanes/local/include -L./ -I./ 

#newTest: newTest.cpp ArithmeticHeader.cpp Stream.cpp
#	icpc $(CFLAGS) $(LFQHULL) -mkl -qopenmp -c -std=c++11

#foo: ArithmeticHeader.cpp Stream.cpp
#	icpc $(CFLAGS) -o newTest newTest.cpp -mkl -qopenmp -pg -std=c++11


