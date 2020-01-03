#pragma once
#include <iostream>
#include "ArithmeticHeader.h"
using namespace std;

class Polyhedron{
	double* pointArray;
	Triangle *triangleArray;
	int tricount;
	int extcount;
	int *exteriorArray;
	int iteration;

public:
	Polyhedron(string input,string output);
	~Polyhedron();
	void sampleLargeTriangles();
	double distance(double x1, double y1, double z1, double x2, double y2, double z2);
	void setTriCount(string input,string output);
	bool isInterior(Triangle* tri);
	void copyToExterior();
	bool pointInTriangle(Triangle* tri,double x,double y,double z);
	void parseToOutput(string input, string output);
	double*  readOutputToArrays(string output);
	bool checkDistanceCutoff(int index, double *b, Triangle* tri);
	void combineDocs();
};
