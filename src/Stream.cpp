#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "StreamHeader.h"
#include <streambuf>
#include <istream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <bits/stdc++.h>
#include "cySampleElim.h"
#include "cyPointCloud.h"
#include "cyHeap.h"
#include "cyPoint.h"
//#include "ArithmeticHeader.h"
#include <cstdlib>
#include <ctime>
#include <stdexcept>
/*
//includes related to C++
#include "libqhull_r/libqhull_r.h"
//#include "libqhullcpp/Qhull.h"

#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullSet.h"
#include "libqhullcpp/QhullVertexSet.h"
#include <vector>
*/
using std::cerr;
using std::cin;
using std::cout;
using std::endl;

//using namespace std;

Polyhedron::Polyhedron(string input, string output)
{
	int pointct = 0;
	parseToOutput(input,output);
	string line;
	ifstream file("../data/" + output);
	while(getline(file,line)){
		pointct++;
	}
	double *points = new double[3*pointct-3];
	points = readOutputToArrays(output);
	pointArray = points;
        //coordT *pointsT;
	//orgQhull::Qhull qhull;
        //memcpy(pointsT,points,(3*pointct-3));
	Triangle *triArray = new Triangle[(pointct-1)/3];
	setTriCount(input,output);
	for(int i = 0; i < (3*pointct-3);i+=9){
		double *aa = new double[3] {pointArray[i+0],pointArray[i+1],pointArray[i+2]};
		double *bb = new double[3] {pointArray[i+3],pointArray[i+4],pointArray[i+5]};
		double *cc = new double[3] {pointArray[i+6],pointArray[i+7],pointArray[i+8]};
		Triangle *bee = new Triangle(aa,bb,cc);
		triArray[i/9] = *bee;
	}
	triangleArray = triArray;
	copyToExterior();
	iteration = 0;
	sampleLargeTriangles();
	//delaunay, creates inputSol
	//parseToOutput(inputSol,outputSol)
	//reset points array and tri array
	//copyToExterior()
	//voronoiLargeTriangles
	//makeGaussianDoc()
	cout << "MADE A POLYHEDRON" << endl;
}

Polyhedron::~Polyhedron()
{
	delete[] pointArray;
	delete[] exteriorArray;
	delete[] triangleArray;
	cout<< "DELETED A POLYHEDRON" << endl;
}

void Polyhedron::setTriCount(string input, string output){
	int pointct = 0;
        parseToOutput(input,output);
        string line;
        ifstream file("../data/" + output);
        while(getline(file,line)){
                pointct++;
        }
	tricount = ((pointct-1)/3);
}

bool Polyhedron::pointInTriangle(Triangle* tri, double x, double y ,double z){
	bool boolz = false;
	double *points;
	points = new double[9];
	tri->getPoints2(points);
	double pt[3] = {x,y,z};
	double pp1[3] = {points[0],points[1],points[2]};
	double pp2[3] = {points[3],points[4],points[5]};
	double pp3[3] = {points[6],points[7],points[8]};
	Triangle *T1 = new Triangle(pp1,pp2,pt);
	Triangle *T2 = new Triangle(pp2,pp3,pt);
	Triangle *T3 = new Triangle(pp1,pp3,pt);
	double A1 = T1->getArea();
	double A2 = T2->getArea();
	double A3 = T3->getArea();
	double myArea = tri->getArea();
	double areasum = A1 + A2 + A3;
	double areadiff = areasum - myArea;
	if(abs(areadiff) < 0.005){
		boolz = true;
	}
	delete[] points;
	return boolz;
}

bool Polyhedron::isInterior(Triangle *tri){
        bool flag = true;
        int upNormalCt=0;
        int downNormalCt=0;
	double *triangleInTriList;
	triangleInTriList = new double[9];
	double *points;
	points = new double[9];
	int i;
	//#pragma omp parallel for private(i)
	for(int i = 0; i < tricount; i++){
		triangleArray[i].getPoints2(triangleInTriList);
		tri->getPoints2(points);
	
                double planePoint[3] = {triangleInTriList[0],triangleInTriList[1],triangleInTriList[2]};
                double planeNormal[3] ={triangleArray[i].getNormal()[0],triangleArray[i].getNormal()[1],triangleArray[i].getNormal()[2]};
                double linePoint[3] = {tri->getCenter()[0],tri->getCenter()[1],tri->getCenter()[2]};

                double V1[3] = {points[6]-points[0],points[7]-points[1],points[8]-points[2]};
                double V2[3] = {points[3]-points[0],points[4]-points[1],points[5]-points[2]};
                double lineDir[3] = {V1[1]*V2[2]-V1[2]*V2[1],-(V1[0]*V2[2]-V1[2]*V2[0]),V1[0]*V2[1]-V1[1]*V2[0]};

                double dMag = sqrt(lineDir[0]*lineDir[0] + lineDir[1]*lineDir[1] + lineDir[2]*lineDir[2]);
                double lineDirNormalize[3] = {lineDir[0]/dMag,lineDir[1]/dMag,lineDir[2]/dMag};

                double tNum = (planeNormal[0]*planePoint[0] + planeNormal[1]*planePoint[1] + planeNormal[2]*planePoint[2])
                - (planeNormal[0]*linePoint[0] + planeNormal[1]*linePoint[1] + planeNormal[2]*linePoint[2]);
                double tDenom = planeNormal[0]*lineDirNormalize[0] + planeNormal[1]*lineDirNormalize[1] + planeNormal[2]*lineDirNormalize[2];
                double t = tNum/tDenom;
                double projectedPoint[3] = {linePoint[0] + (lineDirNormalize[0]*t) ,linePoint[1] +(lineDirNormalize[1]*t) ,linePoint[2] +(lineDirNormalize[2]*t)};

                double p1[3] = {triangleInTriList[0],triangleInTriList[1],triangleInTriList[2]};
                double p2[3] = {triangleInTriList[3],triangleInTriList[4],triangleInTriList[5]};
                double p3[3] = {triangleInTriList[6],triangleInTriList[7],triangleInTriList[8]};
                Triangle *tria = new Triangle(p1,p2,p3);
                        if(pointInTriangle(tria,projectedPoint[0],projectedPoint[1],projectedPoint[2]) == true){
                        	double C2[3] = {triangleArray[i].getCenter()[0],triangleArray[i].getCenter()[1],triangleArray[i].getCenter()[2]};
                        	double TVector[3] = {C2[0]-linePoint[0],C2[1]-linePoint[1],C2[2]-linePoint[2]};
                        	double TVectorMag = sqrt(TVector[0]*TVector[0]+TVector[1]*TVector[1]+TVector[2]*TVector[2]);
                        	double UnitT[3] = {TVector[0]/TVectorMag,TVector[1]/TVectorMag,TVector[2]/TVectorMag};
                        
                        	double dot = UnitT[0]*lineDirNormalize[0] + UnitT[1]*lineDirNormalize[1]+UnitT[2]*lineDirNormalize[2];
                                if(dot >= 0){
                                	upNormalCt++;
                        	}
                        	else{
                                	downNormalCt++;
                        	}
                        }
        }
        upNormalCt--;
        downNormalCt--;
	delete[] triangleInTriList;
	delete[] points;
        if(upNormalCt*downNormalCt == 0){
                flag = false;
        }
                return flag;
}
             


void Polyhedron::copyToExterior(){
    int *oneHot = new int[tricount];
    int exteriorCount = 0;
    int i;
    double *pts;
    pts = new double[9];
    //#pragma omp parallel for private(i)
    for(int i = 0; i < (tricount); i++){
		triangleArray[i].getPoints2(pts);
        	double p1[3] = {pts[0],pts[1],pts[2]};
                double p2[3] = {pts[3],pts[4],pts[5]};
                double p3[3] = {pts[6],pts[7],pts[8]};
                Triangle *tria = new Triangle(p1,p2,p3);
            	oneHot[i] = !isInterior(tria);
    }
    for(int j = 0; j < (tricount); j++){
   	exteriorCount += oneHot[j];
    }
    exteriorArray = oneHot;
    delete[] pts;
}


bool Polyhedron::checkDistanceCutoff(int index, double *b, Triangle *tri){
        bool flag = true;
        int upNormalCt=0;
        int downNormalCt=0;
        double *triangleInTriList;
        triangleInTriList = new double[9];
        double *points;
        points = new double[9];
        int i;
        //#pragma omp parallel for private(i)
        for(int i = 0; i < tricount; i++){
		if(exteriorArray[i] == 1){
                	triangleArray[i].getPoints2(triangleInTriList);
			tri->getPoints2(points);		

                double planePoint[3] = {triangleInTriList[0],triangleInTriList[1],triangleInTriList[2]};
                double planeNormal[3] ={triangleArray[i].getNormal()[0],triangleArray[i].getNormal()[1],triangleArray[i].getNormal()[2]};
                double linePoint[3] = {b[0],b[1],b[2]};

                double V1[3] = {points[6]-points[0],points[7]-points[1],points[8]-points[2]};
                double V2[3] = {points[3]-points[0],points[4]-points[1],points[5]-points[2]};
                double lineDir[3] = {V1[1]*V2[2]-V1[2]*V2[1],-(V1[0]*V2[2]-V1[2]*V2[0]),V1[0]*V2[1]-V1[1]*V2[0]};

                double dMag = sqrt(lineDir[0]*lineDir[0] + lineDir[1]*lineDir[1] + lineDir[2]*lineDir[2]);
                double lineDirNormalize[3] = {lineDir[0]/dMag,lineDir[1]/dMag,lineDir[2]/dMag};

                double tNum = (planeNormal[0]*planePoint[0] + planeNormal[1]*planePoint[1] + planeNormal[2]*planePoint[2])
                - (planeNormal[0]*linePoint[0] + planeNormal[1]*linePoint[1] + planeNormal[2]*linePoint[2]);
                double tDenom = planeNormal[0]*lineDirNormalize[0] + planeNormal[1]*lineDirNormalize[1] + planeNormal[2]*lineDirNormalize[2];
                double t = tNum/tDenom;
                double projectedPoint[3] = {linePoint[0] + (lineDirNormalize[0]*t) ,linePoint[1] +(lineDirNormalize[1]*t) ,linePoint[2] +(lineDirNormalize[2]*t)};

                double p1[3] = {triangleInTriList[0],triangleInTriList[1],triangleInTriList[2]};
                double p2[3] = {triangleInTriList[3],triangleInTriList[4],triangleInTriList[5]};
                double p3[3] = {triangleInTriList[6],triangleInTriList[7],triangleInTriList[8]};
                Triangle *tria = new Triangle(p1,p2,p3);
                        if(pointInTriangle(tria,projectedPoint[0],projectedPoint[1],projectedPoint[2]) == true){
                                double C2[3] = {triangleArray[i].getCenter()[0],triangleArray[i].getCenter()[1],triangleArray[i].getCenter()[2]};
                                double TVector[3] = {C2[0]-linePoint[0],C2[1]-linePoint[1],C2[2]-linePoint[2]};
                                double TVectorMag = sqrt(TVector[0]*TVector[0]+TVector[1]*TVector[1]+TVector[2]*TVector[2]);
                                double UnitT[3] = {TVector[0]/TVectorMag,TVector[1]/TVectorMag,TVector[2]/TVectorMag};

                                double dot = UnitT[0]*lineDirNormalize[0] + UnitT[1]*lineDirNormalize[1]+UnitT[2]*lineDirNormalize[2];
                                if(dot >= 0){
                                        upNormalCt++;
                                }
                                else{
                                        downNormalCt++;
                                }
                        }
		}
        }
        upNormalCt--;
        downNormalCt--;
        if(upNormalCt*downNormalCt == 0){
                flag = false;
        }
	else{
		exteriorArray[index] = 0;
	}
	
	delete[] triangleInTriList;
	delete[] points;
	return flag;
}


void Polyhedron::sampleLargeTriangles(){
    int bigCount = 0;
    double bigSurfaceArea = 0;
    double *pts;
    pts = new double[9];
    int i;
    int pointCt=0;
    double *check = new double[3] {1.2, 2.2, 3.3};
    //#pragma omp parallel for private(i)
    for(int i = 0;i < tricount; i++){
	if(exteriorArray[i] == 1){
		triangleArray[i].getPoints2(pts);
        	double p1[3] = {pts[0],pts[1],pts[2]};
                double p2[3] = {pts[3],pts[4],pts[5]};
                double p3[3] = {pts[6],pts[7],pts[8]};
                Triangle *tri = new Triangle(p1,p2,p3);
        	double area = tri->getArea();
		int cumsum = 0;
        	if(area > 3.5){
			for (int k = 0; k < 20; k++){
				check = tri->PlacePointInTriangle();
				cumsum += checkDistanceCutoff(i,check,tri);
				delete[] check;
			}

			if (cumsum == 0){
				bigSurfaceArea+=area;
				for(int i = 0; i < area;i++){
					pointCt++;
				}
            		bigCount++;
			}
			
        	}
	}
    }
	
    double *tempPt;
    std::vector< cy::Point3f > inputPoints(pointCt);
    int ct = 0;
    
    for(int l = 0;l < tricount; l++){
        if((exteriorArray[l] == 1)){
                triangleArray[l].getPoints2(pts);
                double p1[3] = {pts[0],pts[1],pts[2]};
                double p2[3] = {pts[3],pts[4],pts[5]};
                double p3[3] = {pts[6],pts[7],pts[8]};
                Triangle *tri = new Triangle(p1,p2,p3);
                double area = tri->getArea();
		
                if(area > 3.5){
			for(int j = 0; j < area ; j++){
				tempPt = tri->PlacePointInTriangle();
				if (ct < pointCt){
                        		inputPoints[ct].x = (float)tempPt[0];
                        		inputPoints[ct].y = (float)tempPt[1];
                        		inputPoints[ct].z = (float)tempPt[2];
				}
				ct++;
			}
                }
	
        }
    }

    int numpts = floor(bigSurfaceArea/5);
    cy::WeightedSampleElimination< cy::Point3f, float, 2, int > wse;
    vector< cy::Point3f > outputPoints(numpts);
    float d_max = 1.0;
    wse.Eliminate( inputPoints.data(), inputPoints.size(),outputPoints.data(), outputPoints.size(), d_max, 2);

    if(numpts > 0){
    	FILE* fp;
    	fp = fopen("solutionset.txt","w");
    	for(int i = 0; i < numpts; i++){
         	cout << "O     " << outputPoints[i].x << "     " << outputPoints[i].y << "      " << outputPoints[i].z << endl;
         	fprintf(fp,"       %1.6f %1.6f %1.6f\n",outputPoints[i].x,outputPoints[i].y,outputPoints[i].z);
    	}
    	fclose(fp);
    }

    double smallest = 20.0;    
    double counter = 0;
    for(int n = 0; n < tricount; n++){
	for(int j = 0; j < numpts;j++){
		triangleArray[n].getPoints2(pts);
        	double p1[3] = {pts[0],pts[1],pts[2]};
        	double p2[3] = {pts[3],pts[4],pts[5]};
        	double p3[3] = {pts[6],pts[7],pts[8]};
        	Triangle *tri = new Triangle(p1,p2,p3);
		double area = tri->getArea();
		if (area < smallest){
			smallest = area;
		}
		if((exteriorArray[n] == 1) && 
		(pointInTriangle(tri,outputPoints[j].x,outputPoints[j].y,outputPoints[j].z))){
			exteriorArray[n] = 0;
			counter += area;
		}
	}
    }

    double areaRemaining = bigSurfaceArea - counter;
    cout << "ITERATION: " << iteration << endl;
    cout << "AREA REMAINING: " << areaRemaining << endl;
    cout << "SMALLEST: " << smallest << endl;
    if( areaRemaining > smallest){
	iteration++;
	sampleLargeTriangles();
    }    

    delete[] pts;
    delete[] tempPt;

}

void Polyhedron::makeSolution(){
	
}

//**********************************************************************************************************************************************************************************************************//

double Polyhedron::distance(double x1,double y1,double z1,double x2,double y2,double z2){
        double distance = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
        return distance;
}

void Polyhedron::parseToOutput(string input, string output)
{ 
  ifstream is; 
  is.open("../data/" + input);
  cout << "input is: " << input << endl;

if (is.fail())
{
cout << "\nThe file was not successfully opened for reading."
<< "\nPlease check that the file currently exists.\n\n";
exit(1);
}

  ofstream ofs; 
  ofs.open("../data/" + output, ofstream::out);   
  cout << "output is: " << output << endl;                 
  char c; 
  int line_no = 1; 
  cout << "is the infile open? : " << is.is_open() << endl;
  cout << "is the outfile open? : " << ofs.is_open() << endl;
  while (is.get(c)) 
  {
  	if (c == '\n'){ 
  		line_no++;
	} 
  	if (((line_no-1) % 5 == 2) || ((line_no-1) % 5 == 3) || ((line_no-1) % 5 == 4)){ 
  			ofs << c; 
	}
  }                    
  cout << "GOT TO END" << endl;                                                                                  
  ofs.close();                                                                                                               
  is.close(); 
}

double* Polyhedron::readOutputToArrays(string output){
	std::ifstream inFile;
	inFile.open("../data/" + output);

	std::stringstream strStream;
	strStream << inFile.rdbuf();
	string str = strStream.str();
	
	int j = 1;
	int ct = 0;
	int pointct = 0;
    	string line;
    	ifstream file("../data/" + output);

	//cout << "is the infile open? " << inFile.is_open() << endl;
	//cout << "called readOutputToArrays" << endl;


    	while (getline(file, line)){
        	pointct++;	
 	}
	
	double** triPoints = new double*[pointct];
        for(int i = 0; i < pointct; ++i){
                triPoints[i] = new double[3];
        }

	while( j < str.length()){
    		double x, y, z; 
    		std::size_t offset = j;
    		x = std::stod(&str[offset]); 
    		y = std::stod(&str[offset + 9]); 
    		z = std::stod(&str[offset + 19]);
		double *point = new double[3] {x,y,z}; 
		if (ct < pointct){                
			triPoints[ct] = point;
		} 
		j+=28;
		ct++;
	}

	double *PointList = new double[pointct*3-3];
	int i = 0;
	while( i < (pointct-3)){
	for(int j = 0; j < (3*pointct-3); j+=3){
		PointList[j] = triPoints[i][0];
		PointList[j+1] = triPoints[i][1];
		PointList[j+2] = triPoints[i][2];
	i++;
	}
	}
	for(int i = 0; i < pointct; ++i){
		delete[] triPoints[i];
	}

	return PointList;
}
 

