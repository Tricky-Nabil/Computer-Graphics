#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string>
#include<iostream>

#include<GL/glut.h>
using namespace std;

#define pi (2 * acos(0.0))


double camHeight , camAngle , draw_grid , draw_axis , angle;
int length_of_a_square , total_square;
double speed = 10;
double camRotateSpeed = 2;
double camRotateAngle = 10;

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
		point p1 , p2;
		p1 = axis * *this;
		p2 = axis * p1;

		return *this + p1 * sin(angle * (pi / 180)) + p2 * (1 - cos(angle * (pi / 180)));
	}

	void print(){
		cout<<"x : "<<x<<" , y : "<<y<<" , z : "<<z<<endl;
	}
    
};

point pos(0 , -80 , 40) , u(0 , 0 , 1) , l(0 , 0 , 0);
point d , r ;

void makePoint_d_r(){
	d = (pos - l).normalize();
	r = (u * d).normalize();
}

class ball{

	public:
    const int MANUAL_MODE = 0;
	const int SIMULATION_MODE = 1;

	double radius ,  velocity , boundary;
	int stacks , slices , angle , mode , count_of_rotate;
	point position;
    ball(){
        position = point(0 , 0 , 0);
        radius = 3 ; boundary = 10 ; velocity = .2;
		angle = 45 ; stacks = 20 ; slices = 24 ; count_of_rotate = 0;
		mode = MANUAL_MODE;
    }

    void draw(){
        drawBall(); // draw a ball
        drawArrow(); // draw an arrow
    }

	void makeUpperHemiSphere(int i , int j , point points[][100]){
		glBegin(GL_QUADS);{
					//upper hemisphere
				glVertex3f(points[i][j].x , points[i][j].y , points[i][j].z);
				glVertex3f(points[i][j+1].x , points[i][j+1].y , points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x , points[i+1][j+1].y , points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x , points[i+1][j].y , points[i+1][j].z);
			}glEnd();

	}

	void makeLowerHemiSphere(int i , int j , point points[][100]){
		glBegin(GL_QUADS);{
					//lower hemisphere
				glVertex3f(points[i][j].x , points[i][j].y , -points[i][j].z);
				glVertex3f(points[i][j+1].x , points[i][j+1].y , -points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x , points[i+1][j+1].y , -points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x , points[i+1][j].y , -points[i+1][j].z);
			}glEnd();

	}

	void drawBallHelper(){
		point points[100][100];
		for(int i = 0 ; i <= stacks ; i++){
			double h = radius * sin(((double)i / (double)stacks) * (pi / 2));
			double r = radius * cos(((double)i / (double)stacks) * (pi / 2));
			for(int j = 0 ; j <= slices ; j++){
				points[i][j].x = r * cos(((double)j / (double)slices) * 2 * pi);
				points[i][j].y = r * sin(((double)j / (double)slices) * 2 * pi);
				points[i][j].z = h;
			}
		}
		for(int i = 0 ; i < stacks ; i++){
			bool color = true;
			for(int j = 0 ; j < slices ; j++){
				if(color){
					glColor3f(1,0,0);
				}else{
					glColor3f(0,1,0);
				}
				makeUpperHemiSphere(i , j , points);
				
				if(color){
					glColor3f(0,1,0);
				}else{
					glColor3f(1,0,0);
				}

				makeLowerHemiSphere(i , j , points);

				
				if(j % (slices / 8) == 0){
                    color = (!color);
                } 
			}
		}
    }

    void drawBall(){
        glPushMatrix();
		glTranslatef(position.x , position.y , position.z);
		glRotatef(-count_of_rotate % 360 , sin(angle * pi / 180) , -cos(angle * pi / 180) , 0);
		drawBallHelper();
		glPopMatrix();
    }

    void drawArrow(){
        glColor3f(0 , 0 , 1);
		glPushMatrix();
		glTranslatef(position.x , position.y , position.z);
		glRotatef(angle , 0 , 0 , 1);
		glLineWidth(5);
		glBegin(GL_LINES);{
			glVertex3f(0 , 0 , 0);
			glVertex3f(10 , 0 , 0);
		}glEnd();
		glBegin(GL_TRIANGLES);{
			glVertex3f(10 , 0 , 0);
			glVertex3f(8 , 2 , 0);
			glVertex3f(8 , -2 , 0);
		}glEnd();
		glLineWidth(1);
		glPopMatrix();
    }

    

    void moveForward(){
        char inBoundary = isInBoundary(position.x, position.y);
		
		if(inBoundary != 'r'){
			double vx = cos(angle * pi / 180);
			double vy = sin(angle * pi / 180);

			if(inBoundary == 'x') vx = -vx;
			if(inBoundary == 'y') vy = -vy;

			double tempAngle  = atan(vy/vx) / (pi / 180);

			if(vx <= 0) angle = 180 + tempAngle;
			else angle = tempAngle;
		}
	
        double delX = velocity * cos(angle * pi / 180);
        double delY = velocity * sin(angle * pi / 180);

        point temp(delX , delY , 0);
        position = position + temp;

        count_of_rotate = count_of_rotate + 10;
    }

    void moveBackward(){
		char inBoundary = isInBoundary(position.x, position.y);
		
		if(inBoundary != 'r'){
			double vx = cos(angle * pi / 180);
			double vy = sin(angle * pi / 180);

			if(inBoundary == 'x') vx = -vx;
			if(inBoundary == 'y') vy = -vy;

			double tempAngle  = atan(vy/vx) / (pi / 180);

			if(vx <= 0) angle = 180 + tempAngle;
			else angle = tempAngle;
		}

        double delX = (-1.0) * velocity * cos(angle * pi / 180);
        double delY = (-1.0) * velocity * sin(angle * pi / 180);

        point temp(delX , delY , 0);
        position = position + temp;

        count_of_rotate = count_of_rotate - 10;
    }

    

	char isInBoundary(double x , double y , double z = 0){
		if(x + radius >= boundary || x - radius <= -boundary) return 'x';
		else if(y + radius >= boundary || y - radius <= -boundary) return 'y';
		else return 'r';
    }
    
};

ball myBall;

void draw_boundary(double a , double height = 4){

	glBegin(GL_QUADS);{
		glVertex3f(a , a , 0);
		glVertex3f(a , -a , 0);
		glVertex3f(a , -a , height);
		glVertex3f(a , a , height);
	}glEnd();

	glBegin(GL_QUADS);{
		glVertex3f(-a , -a , 0);
		glVertex3f(-a , a , 0);
		glVertex3f(-a , a , height);
		glVertex3f(-a , -a ,height);
	}glEnd();

	glBegin(GL_QUADS);{
		glVertex3f(a , -a , 0);
		glVertex3f(-a , -a , 0);
		glVertex3f(-a , -a , height);
		glVertex3f(a , -a , height);
	}glEnd();

	glBegin(GL_QUADS);{
		glVertex3f(a , a , 0);
		glVertex3f(-a , a , 0);
		glVertex3f(-a , a , height);
		glVertex3f(a , a , height);
	}glEnd();

}

void square_one_side(int i , int length_of_a_square){
	glRotatef(45 , 0 , 0 , 1);
	glTranslatef(i * 2 * length_of_a_square , i * 2 * length_of_a_square , 0);
	glBegin(GL_QUADS);{
		glVertex3f(length_of_a_square , length_of_a_square , 0);
		glVertex3f(length_of_a_square , -length_of_a_square , 0);
		glVertex3f(-length_of_a_square , -length_of_a_square , 0);
		glVertex3f(-length_of_a_square , length_of_a_square , 0);
	}glEnd();
}

void square_other_side(int i , int length_of_a_square){
	glRotatef(45 , 0 , 0 , 1);
	glTranslatef((-i) * 2 * length_of_a_square , (-i) * 2 * length_of_a_square , 0);
	glBegin(GL_QUADS);{
		glVertex3f(length_of_a_square , length_of_a_square , 0);
		glVertex3f(length_of_a_square , -length_of_a_square , 0);
		glVertex3f(-length_of_a_square , -length_of_a_square , 0);
		glVertex3f(-length_of_a_square , length_of_a_square , 0);
	}glEnd();
}

void squares_in_lines(){
	for(int i = 0 ; i < total_square ; i++){		
		glPushMatrix();
		square_one_side(i , length_of_a_square);
		glPopMatrix();   

		glPushMatrix();
		square_other_side(i , length_of_a_square);
		glPopMatrix();   
	}
}

void drawChessBoardHelper1(int i){
	glPushMatrix();
	glTranslatef(sqrt(2) * length_of_a_square + i * 2 * sqrt(2) * length_of_a_square , 0 , 0);
	glColor3f(1.0 , 1.0 , 1.0);
	squares_in_lines();
	glPopMatrix();
}

void drawChessBoardHelper2(int i){
	glPushMatrix();
	glTranslatef((-1) * sqrt(2) * length_of_a_square - i * 2 * sqrt(2) * length_of_a_square , 0 , 0);
	glColor3f(1.0 , 1.0 , 1.0);
	squares_in_lines();
	glPopMatrix();
}

void drawChessBoard(){
	for(int i = 0 ; i < total_square ; i++){
		drawChessBoardHelper1(i);
		drawChessBoardHelper2(i);	
	}
}
void cameraView(){
	gluLookAt(pos.x , pos.y , pos.z,
        l.x , l.y , l.z,
        u.x , u.y , u.z);
}


void keyboardListener(unsigned char key, int x,int y){

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
			// cout<<"Before rotate"<<endl;
			// u.print();
			u = u.rotate(d , camRotateAngle);
			makePoint_d_r();
			// cout<<"After rotate: "<<endl;
			// u.print();
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
		
		case 'j': 
		case 'J':
			myBall.angle = (myBall.angle + 20) % 360;
			break;

		case 'l': 
		case 'L':
			myBall.angle = (myBall.angle - 20 + 360) % 360;
			break;

		case 'i': 
		case 'I':	
			if(myBall.mode == 0){
				myBall.moveForward();
			}
			break;

		case 'k': 
		case 'K':
			if(myBall.mode == 0){
				myBall.moveBackward();
			}
			break;

		case ' ':
			if(myBall.mode == 1) myBall.mode = 0;
			else myBall.mode = 1;
			break;
		
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
				pos = pos + d * speed;
				break;
		case GLUT_KEY_UP:		// up arrow key
				pos = pos - d * speed;
				break;

		case GLUT_KEY_RIGHT:
				pos = pos - r * speed;
				break;
		case GLUT_KEY_LEFT:
				pos = pos + r * speed;
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
		// case GLUT_LEFT_BUTTON:
		// 	if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				
		// 	}
		// 	break;

		// default:
		// 	break;
	}
}

void display(){
    //clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0 , 0 , 0 , 0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    cameraView();
	glMatrixMode(GL_MODELVIEW);

	glPushMatrix();
	glTranslatef(0 , 0 , -3);
	drawChessBoard();
	glPopMatrix();

	// square on grid
	glColor3f(1 , 0 , 0);
	glPushMatrix();
	glTranslatef(0 , 0 , -3);
	glRotatef(45 , 0 , 0 , 1);
	draw_boundary(6 * length_of_a_square);
	glPopMatrix();


	// ball
	myBall.boundary = 6 * length_of_a_square;
	glPushMatrix();
	glRotatef(45, 0, 0, 1);
	myBall.draw();
	glPopMatrix();

	glutSwapBuffers();
}
void animate(){
    if(myBall.mode == 1)
		myBall.moveForward();

	angle = angle + 0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
    camHeight=150.0;camAngle=1.0;draw_grid=0;draw_axis=1;angle=0;

	length_of_a_square = 5;
	total_square = 20;
	makePoint_d_r();

	glClearColor(0 , 0 , 0 , 0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(80 ,	1 ,	1 ,	1000.0);
	
}

int main(int argc, char **argv){
    glutInit(&argc, argv);
    glutInitWindowPosition(100 , 100);
    glutInitWindowSize(700 , 700);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
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