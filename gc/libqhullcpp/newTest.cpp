 /* Includes */
 #include <iostream>
 #include "mkl.h"
 #include <cmath>
 #include "ArithmeticHeader.h"
 #include "ArithmeticHeader.cpp"
 #include <ctime> 
 #include "StreamHeader.h" 
 #include "Stream.cpp"
 //#include "libqhull.h"

#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/Qhull.h"

#include <cstdio>   /* for printf() of help message */
#include <ostream>
#include <stdexcept>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;

using namespace orgQhull;

using orgQhull::Qhull;
using orgQhull::QhullError;  
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullQh;
using orgQhull::RboxPoints;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexSet; 

 using namespace std;
 /* Main program */
 int main ()
 {
	//aaa
	srand(unsigned(time(NULL)));
	Polyhedron *poly = new Polyhedron("Test2.txt","output2.txt");
	delete poly;
	/*
	int ct2 = 0;
	while(ct2<100){
		tri->PlacePointInTriangle();
		ct2++;
	}
	delete tri;
	*/
	/*
	cout << "*******************" << endl;
	//trying to make another triangle
	double *p1A = new double[3] {4.253227,-3.090114,-8.506542};
	double *p2A = new double[3] {-2.763880,-8.506493,-4.472198};
	double *p3A = new double[3] {7.236073,-5.257253,-4.472195};
	
	Triangle *trial2 = new Triangle(p1A,p2A,p3A);
	int ct = 0;
	while(ct<100){
		trial2->PlacePointInTriangle();
		ct++;
	}
	delete trial2;
	cout << "********************" << endl;
	*/
	return 0;
 }
