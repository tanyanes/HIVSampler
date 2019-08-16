#pragma once
#include <iostream>
using namespace std;

class Triangle{
	double *point1;
	double *point2;
	double *point3;
	double *transmatrix;
	double *center;
	double *transvector;
	double *normal;
	double area;

public:
	Triangle();
	Triangle(double*,double*,double*);
	~Triangle();
	double *getNormal();
	double getDistance(double x1, double y1, double z1, double x2, double y2, double z2);
	double *getTransMatrix();
	double getArea();
	double* getCenter();
	void getPoints2(double*);
	double* getPoints();
	void setCenter();
	void setTranslationVector();
	void setTransMatrix();
	double *crossproduct(double*, double*);
	double* PlacePointInTriangle();
};
