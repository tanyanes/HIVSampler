 /* Includes */
 #include <iostream>
 #include "mkl.h"
 #include <cmath>
// #include "ArithmeticHeader.h"
 #include <ctime> 
 #include "StreamHeader.h" 
// #include "libqhull.h"
#include <ostream>

 using namespace std;
 int main ()
 {
	srand(unsigned(time(NULL)));
	cout << "COMPILATION WORKED!!!" << endl;
	cout << "aaaaaa" << endl;
	/*
	cout << "*******************" << endl;
	//trying to make a triangle, fill with 100 points
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
	Polyhedron *poly = new Polyhedron("Test2.txt","output2.txt");
	delete poly;
	return 0;
 }
