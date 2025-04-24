#include<iostream>
#include<algorithm>
#include<cmath>
#include<vector>
#include<array>
#include <GL/glut.h>

using namespace std;

#define pi (2 * acos(0.0))

double camHeight , camAngle;

double speed = 10;
double camRotateSpeed = 2;
double camRotateAngle = 10;


double sideLengthOfOct = 50;
double radiusOfCircle = sideLengthOfOct / sqrt(3);

double radius = 0;
double currentSide = 50;
int counter = 0;
unsigned int slices = 24;
int subdivisions = 6;

int rotation = 0;

class point{
    public:
	double x , y , z;
    
    point() = default;
	point(double a , double b , double c){
        x = a ; y = b ; z = c;
	}

    point operator+(point const& obj){
        point p;
        p.x = x + obj.x;
        p.y = y + obj.y;
        p.z = z + obj.z;
        return p;
    }

    point operator-(point const& obj){
        point p;
        p.x = x - obj.x;
        p.y = y - obj.y;
        p.z = z - obj.z;
        return p;
    }

    point operator*(double a){
        point p;
        p.x = x * a;
        p.y = y * a;
        p.z = z * a;
        return p;
    }

	point operator*(point const& obj){
        point p;
        p.x = y * obj.z - z * obj.y;
        p.y = z * obj.x - x * obj.z;
        p.z = x * obj.y - y * obj.x;
        return p;

    }

    point operator/(int a){
        point p;
        p.x = x / a;
        p.y = y / a;
        p.z = z / a;
        return p;
    }

	double dotProduct(point p){
		return x * p.x + y * p.y + z * p.z;
	}

	double sqrt_of_dot_product(){
		return sqrt(dotProduct(*this));
	}

	point normalize(){
		point p;
		double temp = sqrt_of_dot_product();
		p.x = x / temp;
		p.y = y / temp;
		p.z = z / temp;
		return p;
	}

	point rotate(point axis , double angle){
		point p1 = axis * *this;
		point p2 = axis * p1;

		point temp;
		temp = *this + p1 * sin(angle * (pi / 180)) + p2 * (1 - cos(angle * (pi / 180)));
		return temp;
	}


    
};

point pos(60 , -60 , 60) , u(0 , 0 , 1) , l(0 , 0 , 0);
point d , r ;

void makePoint_d_r(){
	d = (pos - l).normalize();
	r = (u * d).normalize();
}

void drawAxes(){
	glColor3f(1.0 , 0.0 , 0.0);
	glBegin(GL_LINES);{
		glVertex3f(100 , 0 , 0);
		glVertex3f(-100 , 0 , 0);
	}glEnd();

	glColor3f(0.0 , 1.0 , 0.0);
	glBegin(GL_LINES);{
		glVertex3f(0 , -100 , 0);
		glVertex3f(0 , 100 , 0);
	}glEnd();

	glColor3f(0.0 , 0.0 , 1.0);
	glBegin(GL_LINES);{
		glVertex3f(0,0 , 100);
		glVertex3f(0,0 ,-100);
	}glEnd();
}

void drawUpperSurafaceOfOct(int i){
    glPushMatrix();
    glRotatef(i * 90 , 0 , 0 , 1);
	glTranslatef(radius / sqrt(3) , radius / sqrt(3) , radius / sqrt(3));
    glBegin(GL_TRIANGLES);{
        glVertex3f(currentSide , 0 , 0);
        glVertex3f(0 , currentSide , 0);
        glVertex3f(0 , 0 , currentSide);
    }glEnd();
	glPopMatrix();
}

void drawLowerSurfaceOfOct(int i){
    glPushMatrix();
    glRotatef(i * 90 , 0 , 0 , 1);
	glTranslatef(radius / sqrt(3) , radius / sqrt(3) , (-1) * radius / sqrt(3));
    glBegin(GL_TRIANGLES);{
        glVertex3f(currentSide , 0 , 0);
        glVertex3f(0 , currentSide , 0);
        glVertex3f(0 , 0 , -currentSide);
    }glEnd();
    glPopMatrix();
}

void drawOctahedron(){
    glColor3f(1 , 0 , 0);

    // generating 4 surfaces
    for(int i = 0 ; i < 4 ; i++){
        //opposite surface same color
        if(i%2){
            glColor3f(1 , 0 , 1);
        }else{
            glColor3f(0 , 1 , 1);
        }
        drawUpperSurafaceOfOct(i);

        if(i%2){
            glColor3f(0 , 1 , 1);
        }else{
            glColor3f(1 , 0 , 1);
        }
        drawLowerSurfaceOfOct(i);
		
    }
}

vector<point> makePoint(){
	vector<point> points(slices+1);

    for(int i = 0 ; i <= slices ; i++){
        double angle = -45 +  (double)i / slices * 90 ;
        point temp(radius * cos(angle * pi / 180) , radius * sin(angle * pi / 180) , 0);
        points[i] = temp;
    }
	return points;
}

void drawCylinderHelper2(double height){
    if(height <= 0) return ;

	auto points = makePoint();
	
	for(int i = 0 ; i < slices ; i++){
		glBegin(GL_QUADS);{
        	glVertex3f(points[i].x , points[i].y , -height);
        	glVertex3f(points[i].x , points[i].y , height);
        	glVertex3f(points[i+1].x , points[i+1].y , height);
        	glVertex3f(points[i+1].x , points[i+1].y , -height);
    	}glEnd();
	}

}


void drawCylinderHelper(){
    double height = currentSide * sqrt(2) ;
    double centerOfCylinder = currentSide / 2;

    if(height < 0) return ;

    glColor3f(1 , 1 , 0);

	for(int i = 0 ; i < 4 ; i++){
        glPushMatrix();
        glRotatef(i * 90 , 0 , 1 , 0);
        glTranslatef(centerOfCylinder , 0 , (-1) * centerOfCylinder);
        glRotatef(45 , 0 , 1 , 0);
        drawCylinderHelper2(height/2);
        glPopMatrix();
    }
    
}

void drawCylinder(){
    // glColor3f(1,0,0);
    drawCylinderHelper();
    
    // glColor3f(0,1,0);
    glPushMatrix();
    glRotatef(90 , 1 , 0 , 0);
    drawCylinderHelper();
    glPopMatrix();

    // glColor3f(0,0,1);
    glPushMatrix();
    glRotatef(90 , 0 , 0 , 1);
    drawCylinderHelper();
    glPopMatrix();
}


vector<vector<point>> buildUnitPositiveX(){

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = 1 << subdivisions + 1;

    vector<vector<point>> vertices(pointsPerRow, vector<point>(pointsPerRow));
    point n1;        // normal of longitudinal plane rotating along Y-axis
    point n2;        // normal of latitudinal plane rotating along Z-axis
    point v;         // direction vector intersecting 2 planes, n1 x n2
    float a1;           // longitudinal angle along Y-axis
    float a2;           // latitudinal angle along Z-axis

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for(unsigned int i = 0 ; i < pointsPerRow ; i++){
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2

		// Latitudal angle and normal of latitudal plan alnog Z axis
        float a2 = (pi / 180) * (45.0f - 90.0f * i / (pointsPerRow - 1)); 
        point n2 = point(-sin(a2) , cos(a2) , 0);

        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for(unsigned int j = 0 ; j < pointsPerRow ; ++j){        // normal for latitudinal plane

            // normal for longitudinal plane
            // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
            // therefore, it is rotating (0,0,-1) vector by longitude angle a1

			// Longitudinal angle and normal of latitudal plan alnog Y axis
            float a1 = (pi / 180) * (-45.0f + 90.0f * j / (pointsPerRow - 1));
            point n1 = point(-sin(a1) , 0 , -cos(a1));

            // find direction vector of intersected line, n1 x n2
            v = (n1 * n2).normalize();

            // add a vertex into array
            vertices[i][j] = v;
        }
    }

    return vertices;
}

void drawSphereHelper(){
	auto vertices = buildUnitPositiveX();

	int pointsPerRow = vertices.size();

	for(int i = 0 ; i < pointsPerRow - 1 ; i++){
		for(int j = 0 ; j < pointsPerRow - 1 ; j++){
            glBegin(GL_QUADS);{
				glVertex3f(vertices[i][j].x , vertices[i][j].y , vertices[i][j].z);
				glVertex3f(vertices[i][j+1].x , vertices[i][j+1].y , vertices[i][j+1].z);
				glVertex3f(vertices[i+1][j+1].x , vertices[i+1][j+1].y , vertices[i+1][j+1].z);
				glVertex3f(vertices[i+1][j].x , vertices[i+1][j].y , vertices[i+1][j].z);
			}glEnd();	
		}
	}
}

void drawFirstFace(){
    glPushMatrix();
	glTranslatef(currentSide , 0 , 0);
	glScalef(radius , radius , radius);
	drawSphereHelper();
	glPopMatrix();
}

void drawSecondFace(){
    glPushMatrix();
	glTranslatef(-currentSide , 0 , 0);
	glScalef(radius , radius , radius);
	glRotatef(180, 0, 0, 1);
	drawSphereHelper();
	glPopMatrix();
}

// draw two faces in one axis
void drawSphereAxis(){
    drawFirstFace();
    drawSecondFace();
}

void drawSphere(){

    glColor3f(1 , 0 , 0);
    drawSphereAxis();

    glRotatef(90 , 0 , 1 , 0);
    glColor3f(0 , 0 , 1);
    drawSphereAxis();

    glRotatef(90 , 0 , 0 , 1);
    glColor3f(0 , 1 , 0);
    drawSphereAxis();

}

void keyboardListener(unsigned char key, int x,int y){
	radius = counter / 100.0 * radiusOfCircle;
	currentSide = sideLengthOfOct - counter / 100.0 * sideLengthOfOct;
	
    switch(key){
        case '1':
			l = l - r * camRotateSpeed;
			makePoint_d_r();
			break;
		case '2':
			l = l + r * camRotateSpeed;
			makePoint_d_r();
			break;
		case '3':
			l = l + u * camRotateSpeed;
			makePoint_d_r();
			break;
		
		case '4':
			l = l - u * camRotateSpeed;
			makePoint_d_r();
			break;

		case '5':
			u = u.rotate(d , camRotateAngle);
			makePoint_d_r();
			break;

		case '6':
			u = u.rotate(d , -camRotateAngle);
			makePoint_d_r();
			break;

		case 'w':
			pos = pos + u * speed;
			break;
		case 's':
			pos = pos - u * speed;
			break;

        case ',': 
        case '<':
            counter = min(counter + 5 , 100);
            break;

        case '.': 
        case '>':
            counter = max(counter - 5 , 0);
            break;
		
		case 'a': 
        case 'A':
			rotation = (rotation - 2 + 360) % 360;
			break;

		case 'd': 
        case 'D':
			rotation = (rotation + 2) % 360;
			break;
		
		default:
			break;
	}

	
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			pos = pos + d * speed;
            camHeight -= 3.0;
			break;
		case GLUT_KEY_UP:		// up arrow key
			pos = pos - d * speed;
            camHeight += 3.0;
			break;

		case GLUT_KEY_RIGHT:
            pos = pos - r * speed;  
			camAngle += 0.03;
			break;
		case GLUT_KEY_LEFT:
            pos = pos + r * speed;
			camAngle -= 0.03;
			break;

		case GLUT_KEY_PAGE_UP:
            pos = pos + u * speed;
			l = l + u * speed;
			break;
		case GLUT_KEY_PAGE_DOWN:
            pos = pos - u * speed;
			l = l - u * speed;
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}

void cameraView(){
	gluLookAt(pos.x , pos.y , pos.z,
        l.x , l.y , l.z,
        u.x , u.y , u.z);
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
	
	// gluLookAt(100*cos(camAngle), 100*sin(camAngle), camHeight,		0,0,0,		0,0,1);
	cameraView();

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);

	//add everything
	drawAxes();
	glPushMatrix();

	glRotatef(rotation , 0 , 0 , 1);

    drawOctahedron();

    drawCylinder();

	drawSphere();

	glPopMatrix();

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	camHeight = 150.0 ; camAngle = 1.0;

	counter = 0;
	radius = (counter / 100.0) * radiusOfCircle;
	currentSide = sideLengthOfOct - (counter / 100.0) * sideLengthOfOct;

    makePoint_d_r();

	//clear the screen
	glClearColor(0 , 0 , 0 , 0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80 ,	1 ,	1 ,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(650 , 650);
	glutInitWindowPosition(100 , 100);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("OpenGL Demo");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}