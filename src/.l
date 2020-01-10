LFQHULL=-L/home/tanyanes/local/lib -lqhull
#LDFLAGS=$(LFQHULL)
CFLAGS=-I/home/tanyanes/local/include -L./ -I./ 

newTest: newTest.cpp ArithmeticHeader.cpp Stream.cpp
	icpc $(CFLAGS) $(LFQHULL) -mkl -qopenmp -c -std=c++11

foo: ArithmeticHeader.cpp Stream.cpp
	icpc $(CFLAGS) -o newTest newTest.cpp -mkl -qopenmp -pg -std=c++11


