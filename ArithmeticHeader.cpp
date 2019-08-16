  #include <cmath>
  #include <iostream>
  #include <stdio.h> 
  #include <ctime>
  #include "ArithmeticHeader.h"
  #include "mkl.h"
  using namespace std;

Triangle::Triangle(){
	double *a = new double[3] {0.0,0.0,0.0};
	double *b = new double[3] {1.0,0.0,0.0};
	double *c = new double[3] {0.5,(sqrt(3)/2),0.0};
	point1 = a;
	point2 = b;
	point3 = c;
	setTransMatrix();
	double *defcenter = new double[3] {0.5,(sqrt(3)/6),0};
	center = defcenter;
	setTranslationVector();	
	double V1[3] = {point3[0]-point1[0],point3[1]-point1[1],point3[2]-point1[2]};
        double V2[3] = {point2[0]-point1[0],point2[1]-point1[1],point2[2]-point1[2]};
        normal = crossproduct(V1,V2);
	area = 0.5*sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
}

Triangle::Triangle(double* p1, double* p2, double* p3){
	point1 = p1;
	point2 = p2;
	point3 = p3;
	setTransMatrix();
	setCenter();
	setTranslationVector();
	double V1[3] = {point3[0]-point1[0],point3[1]-point1[1],point3[2]-point1[2]};
        double V2[3] = {point2[0]-point1[0],point2[1]-point1[1],point2[2]-point1[2]};
        normal = crossproduct(V2,V1);
        area = 0.5*sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
	//cout << "NORMAL: {" << normal[0] << " , " << normal[1] << " , " << normal[2] << "}" << endl;
	//cout << "AREA: " << area << endl; 
	//cout << "TRIANGLE MADE!!!" << endl;
}

double* Triangle::getNormal(){
	//cout << point1[0] << "  " << point1[1] << "  " << point1[2] << endl;
	//cout << point2[0] << "  " << point2[1] << "  " << point2[2] << endl;
	//cout << point3[0] << "  " << point3[1] << "  " << point3[2] << endl;
	return normal;
} 

double Triangle::getArea(){
	return area; 
}

double* Triangle::getCenter(){
	return center;
}

void Triangle::getPoints2(double *ptI){
	ptI[0] = point1[0];
	ptI[1] = point1[1];
	ptI[2] = point1[2];
	ptI[3] = point2[0];
	ptI[4] = point2[1];
	ptI[5] = point2[2];
	ptI[6] = point3[0];
	ptI[7] = point3[1];
	ptI[8] = point3[2];
	return;
}

double * Triangle::getPoints(){
	double *points = new double[9] {point1[0],point1[1],point1[2],point2[0],point2[1],point2[2],point3[0],point3[1],point3[2]};
	return points;
}

void Triangle::setCenter(){
	double *cent = new double[3] {(point1[0]+point2[0]+point3[0])/3.0f,(point1[1]+point2[1]+point3[1])/3.0f,(point1[2]+point2[2]+point3[2])/3.0f};
	center = cent;
	//printf("CENTER FOR TARGET TRIANGLE: %1.8f  %1.8f %1.8f \n", center[0],center[1],center[2]);
}

void Triangle::setTranslationVector(){
	double *transCenter = new double[3] {0.0,0.0,0.0};
	double *dpoint = new double[3] {0.3329491,0.3329491,-0.00013};
	int M = 3;
	int N = 1;
	int K = 3;
	double ALPHA = 1.0f;
	double BETA = 0.0f;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,M,N,K,ALPHA,transmatrix,K,dpoint,N,BETA,transCenter,N);
	double *tvector = new double[3] {center[0] - transCenter[0], center[1]-transCenter[1], center[2]-transCenter[2]};
	transvector = tvector;
}

double* Triangle::getTransMatrix(){
	return transmatrix;
}

Triangle::~Triangle(){
	delete[] point1;
	delete[] point2;
	delete[] point3;
	delete[] center;
	delete[] transmatrix;
	delete[] transvector;
	delete[] normal;
}

void Triangle::setTransMatrix(){
	double V1[3] = {point3[0]-point1[0],point3[1]-point1[1],point3[2]-point1[2]};
	double V2[3] = {point2[0]-point1[0],point2[1]-point1[1],point2[2]-point1[2]};
	double *crossV1V2 = crossproduct(V1,V2);
	double *A = new double[9] {V1[0],V2[0],crossV1V2[0],V1[1],V2[1],crossV1V2[1],V1[2],V2[2],crossV1V2[2]};
	delete[] crossV1V2;
	transmatrix = A;
}

double* Triangle::crossproduct(double *a, double *b){
	double *output = new double[3];
	output[0] = a[1]*b[2] - a[2]*b[1];
	output[1] = -(a[0]*b[2] - a[2]*b[0]);
	output[2] = a[0]*b[1] - a[1]*b[0];
	return output;
}

double Triangle::getDistance(double x1,double y1,double z1,double x2,double y2,double z2){
	double distance = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
	return distance;
}

double* Triangle::PlacePointInTriangle(){
	//getting random point out of unit square
	double x = 1.2f;
	double y = 2.2f;
	while((x+y) > 1){
		x = (rand() % 10000000) / 10000000.0f;
		y = (rand() % 10000000) / 10000000.0f;
	}	
	double *point = new double[3] {x,y,0.0};

	int M = 3; //n_rows_in_A
	int N = 1; //n_columns_in_C
	int K = 3; //n_columns_in_A
	double ALPHA = 1.0f;
	double BETA = 0.0f;
	double *C = new double[3] {0.0,0.0,0.0};
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, ALPHA, transmatrix, K, point, N, BETA, C, N);
	double *translated = new double[3] {C[0] + transvector[0] , C[1] + transvector[1] , C[2] + transvector[2]}; 

	double r = 1.0f;
	double para1 = getDistance(translated[0],translated[1],translated[2],point1[0],point1[1],point1[2]);
	double para2 = getDistance(translated[0],translated[1],translated[2],point2[0],point2[1],point2[2]);
	double para3 = getDistance(translated[0],translated[1],translated[2],point3[0],point3[1],point3[2]);
	
	//cout << para1 << "  |  " << para2 << "  |  " << para3 << endl;
	
	if((para1 < r) || (para2 < r) || (para3 < r)){
		PlacePointInTriangle();
	}
	else{
		printf("O               %1.8f                %1.8f            %1.8f\n",translated[0],translated[1],translated[2]);
		return translated;
	}
	delete[] C;
	delete[] translated;
}
