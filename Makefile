#LFQHULL=-L/home/tanyanes/local/lib -lqhull
#LDFLAGS=$(LFQHULL)
CFLAGS=-I/home/tanyanes/local/include 

newTest: newTest.cpp ArithmeticHeader.cpp Stream.cpp
	icpc $(CFLAGS) -o newTest newTest.cpp -mkl -qopenmp -pg -std=c++11
