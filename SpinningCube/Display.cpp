#include "Display.h"

Display::Display() {
	// The default values for the camera are 0, -30, 25
	// for the purpose of this project, but it can be changed
	camPoint[0] = 0;
	camPoint[1] = -40;
	camPoint[2] = 25;
	// Default values for directionVector
	directionVector[0] = 0;
	directionVector[1] = 0;
	directionVector[2] = 0;
	//set tValue to 0
	tValue = 0;
	// Default planeConstant is 20 (plane is at y = -20)
	planeConstant = 25;
	// screen values are ' ' by default
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 40; j++) {
			screen[i][j] = ' ';
		}
	}
}

void Display::setTValue(double newTVal) {
	this->tValue = newTVal;
}

void Display::setCamPoints(double x, double y, double z) {
	camPoint[0] = x;
	camPoint[1] = y;
	camPoint[2] = z;
}

void Display::setPlaneConstant(double newPC) {
	this->planeConstant = newPC;
}

void Display::setDirectionVector(double x, double y, double z) {
	directionVector[0] = x;
	directionVector[1] = y;
	directionVector[2] = z;
}

double Display::getTValue() {
	return this->tValue;
}

double* Display::getCamPointArray() {
	return this->camPoint;
}

double Display::getPlaneConstant() {
	return this->planeConstant;
}

double* Display::getDirectionVector() {
	return this->directionVector;
}

//maps the 3D points to a 2D display
void Display::ThreeDToTwoD(Cube* changeCube) {
	double tempX = 0;
	double tempY = 0;
	double tempZ = 0;

	tempX = 0 - camPoint[0];
	tempY = 0 - camPoint[1];
	tempZ = 0 - camPoint[2];

	setDirectionVector(tempX, tempY, tempZ);

	tValue = -1 * (camPoint[1] + planeConstant) / tempY;

	changeCube->setCenter(camPoint, tValue, directionVector[0], directionVector[2]);


	for (int i = 0; i < 8; i++) {
		tempX = changeCube->getXCoordinates()[i] - camPoint[0];
		tempY = changeCube->getYCoordinates()[i] - camPoint[1];
		tempZ = changeCube->getZCoordinates()[i] - camPoint[2];

		setDirectionVector(tempX, tempY, tempZ);

		tValue = -1 * (camPoint[1] + planeConstant) / tempY;

		changeCube->setCoordinatesTwoDX(i, camPoint[0], tValue, directionVector[0]);
		changeCube->setCoordinatesTwoDZ(i, camPoint[2], tValue, directionVector[2]);
	}
}

void Display::displayCube(Cube cube) {

	sf::RenderWindow window{ sf::VideoMode(1920,1080), "Spinning Cube" };

	double minXValue = 1000, minZValue = 1000, maxXValue = -1000, maxZValue = -1000;


	//find the maximum and minimum x and z values in a rotation
	for (int numRuns = 0; numRuns < 192; numRuns++) {

		if (cube.getRotationCalculation() >= 6.28316) {
			cube.setRotationCalculation(0);
		}

		cube.setCoordinates();
		ThreeDToTwoD(&cube);
		cube.findInvalidPoints(*this);

		//std::cout << newCube.getXCoordinates2D()[1];

		//get the smallest x and z value over the rotation using a bottom point
		if (minXValue > cube.getXCoordinates2D()[1] && cube.getXCoordinates2D()[1] != DBL_MAX) {
			minXValue = cube.getXCoordinates2D()[1];
		}
		if (minZValue > cube.getZCoordinates2D()[1] && cube.getZCoordinates2D()[1] != DBL_MAX) {
			minZValue = cube.getZCoordinates2D()[1];
		}
		//get the largest x and z value over the rotation using a top point
		if (maxXValue < cube.getXCoordinates2D()[0] && cube.getXCoordinates2D()[0] != DBL_MAX) {
			maxXValue = cube.getXCoordinates2D()[0];
		}
		if (maxZValue < cube.getZCoordinates2D()[0] && cube.getZCoordinates2D()[0] != DBL_MAX) {
			maxZValue = cube.getZCoordinates2D()[0];
		}

		cube.setRotationCalculation(cube.getRotationCalculation() + PI / 96);
	}


	for (int numRuns = 0; numRuns < 5000; numRuns++) {

		//double rotationVal = -1;
		//cube.setRotationCalculation(PI / 4);

		if (cube.getRotationCalculation() >= 6.28316) {
			cube.setRotationCalculation(0);
		}

		cube.setCoordinates();
		ThreeDToTwoD(&cube);
		cube.findInvalidPoints(*this);


		//the allParallelograms are made so that indexes 0 and 2 are opposites to each other
		int allParallelograms[6][4] = { {0,1,3,2},{0,1,7,6},{4,5,3,2},{4,5,7,6},{0,2,4,6},{1,3,5,7} };

		cube.scale2DValues(maxZValue, minZValue);

		//this loop makes it so that only valid parallelograms are in the allParallelograms array
		for (int i = 0; i < 6; i++) {
			bool isInvalidPoint = false;
			int invPointIter = 0;
			//for each invalid point, check for it in the array
			while (cube.getInvalidPoints()[invPointIter] != -1) {
				for (int j = 0; j < 4; j++) {
					if (allParallelograms[i][j] == cube.getInvalidPoints()[invPointIter]) {
						isInvalidPoint = true;
					}
				}
				invPointIter++;
			}
			//if the invalid point is contained in the array, disqualify whole sub-array
			if (isInvalidPoint == true) {
				for (int j = 0; j < 4; j++) {
					allParallelograms[i][j] = -1;
				}
			}
		}
		//make 3 4-point shapes
		sf::ConvexShape para1;
		sf::ConvexShape para2;
		sf::ConvexShape para3;

		para1.setPointCount(4);
		para2.setPointCount(4);
		para3.setPointCount(4);

		//run through entire parallelograms array and set their points based on 2D coordinates
		int numVertice = 0;
		for (int i = 0; i < 6; i++) {
			if (allParallelograms[i][0] != -1) {
				//if numVertices is 0, set the points to para1
				if (numVertice == 0) {
					for (int j = 0; j < 4; j++) {
						double zCoord =  -1 * (25 * cube.getZCoordinates2D()[allParallelograms[i][j]]) + 1080;
						para1.setPoint(j, sf::Vector2f(25 * cube.getXCoordinates2D()[allParallelograms[i][j]],
							zCoord));
						//set the color of para1 using getColor functions
						para1.setFillColor(getColor(allParallelograms[i]));
					}
				}
				//if numVertices is 1, set the points to para2
				else if (numVertice == 1) {
					for (int j = 0; j < 4; j++) {
						double zCoord = -1 * (25 * cube.getZCoordinates2D()[allParallelograms[i][j]]) + 1080;
						para2.setPoint(j, sf::Vector2f(25 * cube.getXCoordinates2D()[allParallelograms[i][j]],
							zCoord));
						//set the color of para2 using getColor functions
						para2.setFillColor(getColor(allParallelograms[i]));
					}
				}
				//if numVertices is 2, set the points to para3
				else {
					for (int j = 0; j < 4; j++) {
						double zCoord = -1 * (25 * cube.getZCoordinates2D()[allParallelograms[i][j]]) + 1080;
						para3.setPoint(j, sf::Vector2f(25 * cube.getXCoordinates2D()[allParallelograms[i][j]],
							zCoord));
						//set the color of para3 using getColor functions
						para3.setFillColor(getColor(allParallelograms[i]));
					}
				}
				numVertice++;
			}
		}
		//clear display
		window.clear(sf::Color(0, 0, 0, 255));

		//draw all shapes to the screen
		window.draw(para1);
		window.draw(para2);
		window.draw(para3);

		//display
		window.display();

		//change the angle of the cube by PI / 96
		cube.setRotationCalculation(cube.getRotationCalculation() + PI / 96);
		//approximately 60fps, so wait 17 seconds
		std::this_thread::sleep_for(std::chrono::milliseconds(17));
	}
}

sf::Color Display::getColor(int points[]) {
	int point1 = -1;
	int point2 = -1;
	int point3 = -1;

	point1 = points[0];
	point2 = points[1];
	point3 = points[2];


	if ((point1 % 2) == 0 && (point2 % 2) == 0 && (point3 % 2) == 0) {
		return sf::Color::Green;
	}

	if (point1 % 2 == 1) {
		point1--;
	}
	if (point2 % 2 == 1) {
		point2--;
	}
	if (point3 % 2 == 1) {
		point3--;
	}

	point1 = point1 / 2;
	point2 = point2 / 2;
	point3 = point3 / 2;

	//find the two unique points
	int nonRepPoint1 = -1;
	int nonRepPoint2 = -1;

	nonRepPoint1 = point1;

	if (nonRepPoint1 == point2) {
		nonRepPoint2 = point3;
	}
	else {
		nonRepPoint2 = point2;
	}

	if ((nonRepPoint1 == 0 && nonRepPoint2 == 1) || (nonRepPoint1 == 1 && nonRepPoint2 == 0)) {
		return sf::Color::Red;
	}
	if ((nonRepPoint1 == 2 && nonRepPoint2 == 3) || (nonRepPoint1 == 3 && nonRepPoint2 == 2)) {
		return sf::Color::Red;
	}
	if ((nonRepPoint1 == 0 && nonRepPoint2 == 3) || (nonRepPoint1 == 3 && nonRepPoint2 == 1)) {
		return sf::Color::Blue;
	}
	if ((nonRepPoint1 == 1 && nonRepPoint2 == 2) || (nonRepPoint1 == 2 && nonRepPoint2 == 1)) {
		return sf::Color::Blue;
	}

	return sf::Color::Yellow;
}