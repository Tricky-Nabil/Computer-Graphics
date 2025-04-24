#include<bits/stdc++.h>
#include <GL/glut.h>
#include "1905006_Header.h"
#include "bitmap_image.hpp"

using namespace std;

#define pi (2*acos(0.0))

double cameraHeight , cameraAngle;

int recursionLevel, imageHeight, imageWidth;

bitmap_image image;

//declaring the global variables
vector<Object*> objects;
vector<PointLight*> pointlights;
vector<SpotLight*> spotlights;

//number of images
int imageCount = 1;

//initial values for camera
Point pos(150,70,10);

Point up(0,0,1);
Point rightV(-1 / pow(2,0.5), 1 / pow(2,0.5), 0);
Point look(-1 / pow(2,0.5), -1 / pow(2,0.5), 0);


double speed = 3;
double windowWidth = 500, windowHeight = 500 , viewAngle = 80;

//updating the color value
void updateColorValue(vector<double> &color){
	color[0] = (color[0] > 1) ? 1 : ((color[0] < 0) ? 0 : color[0]);
	color[1] = (color[1] > 1) ? 1 : ((color[1] < 0) ? 0 : color[1]);
	color[2] = (color[2] > 1) ? 1 : ((color[2] < 0) ? 0 : color[2]);
}

/*
-------------------------------------
checking the nearest object
helps to eliminate the far objects
-------------------------------------
*/
int getNearestObjectIndex(Ray &ray, vector<double> &color){
	double tMin = -1;
	int nearestObjectIndex = -1;

	for(int k = 0; k < objects.size(); k++){
		double t = objects[k]->intersect(ray, color, 0);
		if(t > 0){
			if(nearestObjectIndex == -1 || t < tMin)
				tMin = t , nearestObjectIndex = k;
		}	
	}
	return nearestObjectIndex;
}

Point getTopLeft(double planeDistance , double tempWindowHeight , double tempWindowWidth){
	Point tempLook = look * planeDistance;
	Point tempUp = up * tempWindowHeight;
	Point tempRight = rightV * tempWindowWidth;
	Point topLeft = pos + tempLook + tempUp - tempRight;
	return topLeft;
}

vector<double> get_du_dv(){
	vector<double> du_dv(2);
	du_dv[0] = (1.0 * windowWidth) / imageWidth;
	du_dv[1] = (1.0 * windowHeight) / imageHeight;
	return du_dv;
}

// capturing the image
void capture(){

	// initializing the image with black color
	for(int i = 0; i < imageWidth; i++){
		for(int j = 0; j < imageHeight; j++){
			image.set_pixel(i, j, 0, 0, 0);
		}	
	}

    // calculating the plane distance
	double tempWindowHeight = (1.0) * (windowHeight / 2.0);
	double tempWindowWidth = (1.0) * (windowWidth / 2.0);
	double tempAngle = (1.0) * tan((pi / 360) * viewAngle);
	double planeDistance = (1.0) * (tempWindowHeight / tempAngle);

	// calculating the top left point
	Point topLeft = getTopLeft(planeDistance, tempWindowHeight, tempWindowWidth);
	
	// calculating the du and dv
	auto du_dv = get_du_dv();
	double du = du_dv[0];
	double dv = du_dv[1];

	
	Point tempRightDu = rightV * du;
	Point tempUpDv = up * dv;
	topLeft = topLeft + (tempRightDu / 2.0) - (tempUpDv / 2.0);
	// capturing the image
	for(int i = 0; i < imageWidth; i++){
		for(int j = 0; j < imageHeight; j++){
			Point tempRightDuI = tempRightDu * i;
			Point tempUpDvJ = tempUpDv * j;
			Point pixel = topLeft + tempRightDuI - tempUpDvJ;
			Point newPos = pixel - pos;
			Ray ray(pos , newPos);
			vector<double> color(3,0);

			int nearestObjectIndex = getNearestObjectIndex(ray, color);

			if(nearestObjectIndex != -1){
                color[0] = color[1] = color[2] = 0;
				double t = objects[nearestObjectIndex]->intersect(ray, color, 1);
				updateColorValue(color);
				image.set_pixel(i, j, 255 * color[0], 255 * color[1], 255 * color[2]);
			}
		}
	}
	/*
	saving the image
	number count is increasing
	*/
	image.save_image("Output_"+to_string(imageCount)+".bmp");
	imageCount = imageCount + 1;	
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
		case '0':
			// capturing the image
			cout<<"Capturing Image\n";
			capture();
			cout<<"Image Captured\n";
			break;
		case '1':
			// rotating the camera
			rightV = rightV.rotate(up, (pi / 180));
			look = look.rotate(up, (pi / 180));
			break;
		case '2':
			// rotating the camera
			rightV = rightV.rotate(up, -(pi / 180));
			look = look.rotate(up, -(pi / 180));
			break;
		case '3':
			// rotating the camera
			up = up.rotate(rightV, (pi / 180));
			look = look.rotate(rightV, (pi / 180));
			break;
		case '4':
			// rotating the camera
			up = up.rotate(rightV, -(pi / 180));
			look = look.rotate(rightV, -(pi / 180));
			break;
		case '5':
		 	// rotating the camera
			rightV = rightV.rotate(look, (pi / 180));
			up = up.rotate(look, (pi / 180));
			break;
		case '6':
			// rotating the camera
			rightV = rightV.rotate(look, -(pi / 180));
			up = up.rotate(look, -(pi / 180));
			break;
		default:
			// do nothing
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		
			//down arrow key
			pos = pos - look * speed;
			break;
		case GLUT_KEY_UP:		
			// up arrow key
			pos = pos + look * speed;
			break;
		case GLUT_KEY_RIGHT:  
			//right arrow key
			pos = pos + rightV * speed;
			break;
		case GLUT_KEY_LEFT:
			//left arrow key
			pos = pos - rightV * speed;
			break;
		case GLUT_KEY_PAGE_UP:
			//page up
			pos = pos + up * speed;
			break;
		case GLUT_KEY_PAGE_DOWN:
			//page down
			pos = pos - up * speed;
			break;
		default:
			// do nothing
			break;
	}
}

// initializing the camera view
void initCameraView(){
	gluLookAt(pos.getX(), pos.getY(), pos.getZ(), 
			pos.getX() + look.getX(), pos.getY() + look.getY(), pos.getZ() + look.getZ(), 
			up.getX(), up.getY(), up.getZ());
	

}

void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//  gluLookAt(0,0,200,	0,0,0,	0,1,0);
	//gluLookAt(0,100,0,	0,0,0,	0,0,1);

	initCameraView();

	

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects
	for(auto object: objects){
		object->draw();
	}
	//add lights
	for(auto pointlight: pointlights){
		pointlight->draw();
	}
	//add spotlights
	for(auto spotlight: spotlights){
		spotlight->draw();
	}
	/*
	-------------------------------------
	Objects, lights and spotlights are added
	-------------------------------------
	*/

	glutSwapBuffers();
}


void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void loadData(){
	//loading the data from the file
	ifstream in("scene.txt");
	in >> recursionLevel >> imageHeight;

	imageWidth = imageHeight;

	//loading the objects
	int objCount;
	in >> objCount;

	for(int i = 0; i < objCount; i++)
	{
		string objType;
		in >> objType;

		Object *obj;

		if(objType == "sphere"){
			//loading the sphere
            Point center;
            double radius;
            in >> center.x >> center.y >> center.z >> radius;
            vector<double> color(3);
            in >> color[0] >> color[1] >> color[2];
            vector<double> coefficients(4);
            in >> coefficients[0] >> coefficients[1] >> coefficients[2] >> coefficients[3];
            int shine;
            in >> shine;
            obj = new Sphere(center, radius);
            obj->setCoefficients(coefficients);
            obj->setColor(color);
            obj->setShine(shine);
		}
		else if(objType == "triangle"){
			//loading the triangle
            Point a , b , c;
            in >> a.x >> a.y >> a.z >> b.x >> b.y >> b.z >> c.x >> c.y >> c.z;
            vector<double> color(3);
            in >> color[0] >> color[1] >> color[2];
            vector<double> coefficients(4);
            in >> coefficients[0] >> coefficients[1] >> coefficients[2] >> coefficients[3];
            int shine;
            in >> shine;
            obj = new Triangle(a, b, c);
            obj->setCoefficients(coefficients);
            obj->setColor(color);
            obj->setShine(shine);
		}
		else if(objType == "general"){
			//loading the general object
            double a , b , c , d , e , f , g , h , i , j;
            in >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j;
            Point ref;
            in >> ref.x >> ref.y >> ref.z;
            double length , width , height;
            in >> length >> width >> height;
            vector<double> color(3);
            in >> color[0] >> color[1] >> color[2];
            vector<double> coefficients(4);
            in >> coefficients[0] >> coefficients[1] >> coefficients[2] >> coefficients[3];
            int shine;
            in >> shine;
			obj = new General(a, b, c, d, e, f, g, h, i, j);
            obj->setReferencePoint(ref);
            obj->setLength(length);
            obj->setWidth(width);
            obj->setHeight(height);
            obj->setCoefficients(coefficients);
            obj->setColor(color);  
            obj->setShine(shine);

		}
		else{
			cout<<objType<<" is not a valid object type"<<endl;
		}
		objects.push_back(obj);
	}

	int pointlightCount;
	in >> pointlightCount;
	//loading the pointlights
	for(int i = 0; i < pointlightCount; i++){
		PointLight *light = new PointLight();
        Point pos;
        vector<double> color(3);
        in >> pos.x >> pos.y >> pos.z >> color[0] >> color[1] >> color[2];
        light->setPosition(pos);
        light->setColor(color);
		pointlights.push_back(light);
	}

	int spotlightCount;
	in >> spotlightCount;
	//loading the spotlights
	for(int i = 0; i < spotlightCount; i++){
        Point pos;
        vector<double> color(3);
        Point dir;
        double cutoffAngle;
        in >> pos.x >> pos.y >> pos.z >> color[0] >> color[1] >> color[2] >> dir.x >> dir.y >> dir.z >> cutoffAngle;
		SpotLight *spotlight = new SpotLight(pos , color , dir , cutoffAngle);
		spotlights.push_back(spotlight);
	}
	//loading the floor
	Object *floor = new Floor(1000, 20);
	floor->setColor({0.5, 0.5, 0.5});
	floor->setCoefficients({0.4, 0.2, 0.2, 0.2});
	objects.push_back(floor);
	
}


void init(){
	cameraHeight=150.0;
	cameraAngle=1.0;

	loadData();
	image = bitmap_image(imageWidth, imageHeight);

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){

	
	glutInit(&argc,argv);
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(400, 200);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Ray Tracing");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);

	glutMainLoop();		//The main loop of OpenGL

	//deleting the objects
	objects.clear();
	objects.shrink_to_fit();


	pointlights.clear();
	pointlights.shrink_to_fit();

	spotlights.clear();
	spotlights.shrink_to_fit();
	/*
	-------------------------------------
	Objects, lights and spotlights are deleted
	-------------------------------------
	*/

	return 0;
}