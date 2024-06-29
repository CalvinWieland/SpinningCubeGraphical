#pragma once
#include <SFML/Graphics.hpp>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include<chrono> 
#include<thread>
#include <cmath>
#include "Display.h"
#define PI 3.14159
#define INF std::numeric_limits<double>::infinity();

class Display;

class Cube {
private:
	//the x-coordinates and y coordinates are stored in arrays
	double xCoordinates[8];
	double yCoordinates[8];
	double zCoordinates[8];

	double centerOfCubeTwoD[2];

	//the x and y coordinates after 3D projection
	double xCoordinatesTwoD[8];
	double zCoordinatesTwoD[8];

	//make an array which stores which points should not be drawn
	int invalidPoints[4] = {-1};
	
	//rotationSize says how many radians the cube should rotate per rotation
	double rotationSize;
	//rotationCalculation is the current angle of the cube
	double rotationCalculation;

	double cubeLength;

	//x, y, z
	int pointsMap[2][2][2];
	int insNum;

	//tracks which point is being checked for intersections
	bool currentTestPoint[3];

public:
	Cube();

	void rotateCube();
	void findInvalidPoints(Display& camera);

	//Setters
	void setRotationCalculation(double newRotationCalculation);
	void setRotationSize(double newRotationSize);
	void setCoordinates();
	void setCenter(double camPoint[], double tValue, double directionComp, double directionComp2);
	void setCenter(double x, double z);

	//Getters
	double* getXCoordinates();
	double* getYCoordinates();
	double* getZCoordinates();
	double getRotationCalculation();
	double* getXCoordinates2D();
	double* getZCoordinates2D();
	int* getInvalidPoints();
	double* getCenterOfCubeTwoD();

	void setCoordinatesTwoDX(int index, double camPoint, double tValue, double directionComp);
	void setCoordinatesTwoDZ(int index, double camPoint, double tValue, double directionComp);

	//for all perpendicular to z-axis edges, determine if it's on the left, right or middle of axis of roation
	void leftRightOrMiddlePerp(int LRM[], int cubeEdgesPerpendicular[][2]);
	//determines the order of the perpendicluar edges from left to right. this helps determine which symbols to print
	void determineOrder(int perpOrder[], int leftRightorMiddle[]);

	void setCoordinatesTwoDZ(int index, double newValue);
	void setCoordinatesTwoDX(int index, double newValue);

	double getCubeLength();

	void makeParallelograms(std::tuple<bool, bool, bool> para1[], 
		std::tuple<bool, bool, bool> para2[], std::tuple<bool, bool, bool> para3[], 
		bool oppositeToCurPoint[]);

	void scale2DValues(double maxZValue, double minZValue);

	int getPointNum(std::tuple<bool, bool, bool> &testTuple);

	int isCollision(std::tuple<bool, bool, bool> parallelogram[], int testPoint, Display& camera);
};