#include<bits/stdc++.h>
#include "bitmap_image.hpp"
#include <GL/glut.h>

using namespace std;

#define pi (2*acos(0.0))

extern bitmap_image image;

class Point{  
    
    public:
	double x,y,z,t;

    Point() : x(0) , y(0) , z(0) , t(0) {}
	Point(double x , double y , double z) : x(x) , y(y) , z(z) , t(1.0) {}
   
	Point operator+(Point p){
        Point temp;
        temp.x = x + p.x;
        temp.y = y + p.y;
        temp.z = z + p.z;
        return temp;
        
    }
    Point operator-(Point p){
        Point temp;
        temp.x = x - p.x;
        temp.y = y - p.y;
        temp.z = z - p.z;
        return temp;
    }
    Point operator-(){
       Point temp;
        temp.x = -x;   
        temp.y = -y;
        temp.z = -z;
        return temp;
    }
	Point operator*(double val){
        Point temp;
        temp.x = x * val;
        temp.y = y * val;
        temp.z = z * val;
        return temp;
    }
	Point operator/(double b){
        Point temp;
        temp.x = x / b;
        temp.y = y / b;
        temp.z = z / b;
        return temp;
    }
    double dotProduct(Point p){
        double temp;
        temp = x * p.x + y * p.y + z * p.z;
        return temp;
    } 
    Point crossProduct(Point p){
        Point temp;
        temp.x = y * p.z - z * p.y;
        temp.y = z * p.x - x * p.z;
        temp.z = x * p.y - y * p.x;
        return temp;
    } 
    double value(){
        double temp;
        temp = sqrt(dotProduct(*this));
        return temp;
    }
    Point normalize(){
        double len = value();
        x = x / len;
        y = y / len;
        z = z / len;
        return Point(x , y , z);
    }
    Point rotate(Point axis , double angle){
        Point p1;
		p1 = axis.crossProduct(*this);
        return *this * cos(angle) + p1 * sin(angle);
		
	}
    void print(){
		cout<<"x : "<<x<<" , y : "<<y<<" , z : "<<z<<endl;
	}

    double getX(){
        return x;
    }
    double getY(){
        return y;
    }
    double getZ(){
        return z;
    }
    double getT(){
        return t;
    }
    void setX(double x){
        this->x = x;
    }
    void setY(double y){
        this->y = y;
    }
    void setZ(double z){
        this->z = z;
    }
    void setT(double t){
        this->t = t;
    }

    void setValues(double x , double y , double z){
        this->x = x;
        this->y = y;
        this->z = z;
        this->t = 1.0;
    }

};


class PointLight{
    public:
    Point light_pos;
    vector<double> color;

    void draw(){
        glPointSize(5);
        glBegin(GL_POINTS);
        glColor3f(color[0], color[1], color[2]);
        glVertex3f(light_pos.getX(), light_pos.getY(), light_pos.getZ());
        glEnd();
    }

    void setPosition(Point pos){
        this->light_pos.setValues(pos.getX(), pos.getY(), pos.getZ());
    }

    void setColor(vector<double> color){
        this->color = color;
    }

};

class SpotLight{
    public:
    PointLight pointLight;
    Point light_direction;
    double cutoffAngle;

    SpotLight(){
        pointLight.light_pos = Point(0,0,0);
        cutoffAngle = 0;
    }
    SpotLight(Point pos , vector<double> color , Point dir , double cutoffAngle){
        pointLight.light_pos.setValues(pos.getX(), pos.getY(), pos.getZ());
        pointLight.setColor(color);
        this->light_direction.setValues(dir.getX(), dir.getY(), dir.getZ());
        this->cutoffAngle = cutoffAngle;
    }
    void draw(){
        glPointSize(15);
        glBegin(GL_POINTS);
        glColor3f(pointLight.color[0] , pointLight.color[1] , pointLight.color[2]);
        glVertex3f(pointLight.light_pos.getX(), pointLight.light_pos.getY(), pointLight.light_pos.getZ());
        glEnd();
    }

};

class Ray{
    public:
    Point start , dir;
    Ray(Point start, Point dir){
        this->start.setValues(start.getX(), start.getY(), start.getZ());
        Point temp = dir.normalize();
        this->dir.setValues(temp.getX(), temp.getY(), temp.getZ());
    }

};

class Object;

extern vector<PointLight*> pointlights;
extern vector<SpotLight*> spotlights;
extern vector<Object*> objects;
extern int recursionLevel;

class Object {
public:
		Point reference_point;
		double height, width, length;
		vector<double> color;
		vector <double> coefficients; // 0-ambience, 1-diffuse, 2-specular, 3-reflection
		int shine; 
		
		Object(){
            coefficients = vector<double>(4,0);
            color = vector<double>(3,0);
		}

		void setColor(vector<double> color){
            this->color = color;
        }

        virtual vector<double> getColorAt(Point point){
            return this->color;
        }
    	
		void setShine(int shine){
            this->shine = shine;
        }

		void setCoefficients(vector<double> coefficients){
            this->coefficients = coefficients;
        }  
        
        void setReferencePoint(Point reference_point){
            this->reference_point.setValues(reference_point.getX(), reference_point.getY(), reference_point.getZ());
        }

        void setLength(double length){
            this->length = length;
        }

        void setWidth(double width){
            this->width = width;
        }

        void setHeight(double height){
            this->height = height;
        }

        virtual void draw(){}
		virtual double intersectHelper(Ray ray, vector<double> &color, int level) = 0;
        virtual Ray getNormal(Point point, Ray incidentRay) = 0;
        virtual void updateColorWithAmbience(vector<double> &color , vector<double> colorAtIntersection , double ambienceCoefficient){
            color[0] = colorAtIntersection[0] * ambienceCoefficient;
            color[1] = colorAtIntersection[1] * ambienceCoefficient;
            color[2] = colorAtIntersection[2] * ambienceCoefficient;
        }
        virtual void updateDiffuseAndSpecular(vector<double> &color, vector<double> colorAtIntersection, vector<double> lightColor, double val, double phong , double coeff1 , double coeff2){
            double powOfPhong = pow(phong, shine);
            color[0] += lightColor[0] * coeff1 * colorAtIntersection[0] * val;
            color[0] += lightColor[0] * coeff2 * colorAtIntersection[0] * powOfPhong;

            color[1] += lightColor[1] * coeff1 * colorAtIntersection[1] * val;
            color[1] += lightColor[1] * coeff2 * colorAtIntersection[1] * powOfPhong;

            color[2] += lightColor[2] * coeff1 * colorAtIntersection[2] * val;
            color[2] += lightColor[2] * coeff2 * colorAtIntersection[2] * powOfPhong;
        }
        virtual void handleObscured(vector<double> &color, vector<double> colorAtIntersection, Ray lightRay, Ray normal , Point intersectionPoint , Ray ray , int i , bool flag = false){
            Point lightRayDir = lightRay.dir;
            Point normalDir = normal.dir;
            Point rayDir = ray.dir;
            double checkVal = -lightRayDir.dotProduct(normalDir);
            double val;
            if(checkVal < 0){
                val = 0;
            }
            else{
                val = checkVal;
            }

            double temp = 2 * lightRayDir.dotProduct(normalDir);
            Point reflectionDir = lightRayDir - normalDir * temp;
            double checkPhongval = -rayDir.dotProduct(reflectionDir.normalize());
            double phong;
            if(checkPhongval < 0){
                phong = 0;
            }
            else{
                phong = checkPhongval;
            }
            if(flag){
                updateDiffuseAndSpecular(color, colorAtIntersection, pointlights[i]->color, val, phong , coefficients[1] , coefficients[2]);
            }
            else{
                updateDiffuseAndSpecular(color, colorAtIntersection, spotlights[i]->pointLight.color, val, phong , coefficients[1] , coefficients[2]);
            }
        }
        virtual bool checkObscured(Ray lightRay, double t2 , vector<double> &color){
            for(auto obj : objects){
                double t3 = obj->intersectHelper(lightRay, color, 0);
                if(t3 > 0){
                    if(t3 + 1e-5 < t2)
                        return true;
                    else
                        return false;
                }
            }
            return false;
        }
        virtual void calculationsForPointLights(vector<double> &color, vector<double> colorAtIntersection, Point intersectionPoint, Ray ray){
            for(int i = 0 ; i < pointlights.size() ; i++){
                Point lightPosition, lightDirection;
                lightPosition.setValues(pointlights[i]->light_pos.getX(), pointlights[i]->light_pos.getY(), pointlights[i]->light_pos.getZ());
                Point temp = (intersectionPoint - lightPosition).normalize();
                lightDirection.setValues(temp.getX(), temp.getY(), temp.getZ());
               
                Ray lightRay(lightPosition , lightDirection);

                Ray normal = getNormal(intersectionPoint , lightRay);
                
                double t2 = (intersectionPoint - lightPosition).value();
                if(t2 < 1e-5) continue;

                if(!checkObscured(lightRay, t2, color)){
                    handleObscured(color, colorAtIntersection, lightRay, normal , intersectionPoint , ray , i , true); 
                }
                else{
                    /// do nothing
                }
            }
        }
        virtual void calculationsForSpotLights(vector<double> &color, vector<double> colorAtIntersection, Point intersectionPoint, Ray ray){
            for(int i = 0; i < spotlights.size(); i++){
                Point lightPosition,lightDirection;
                lightPosition.setValues(spotlights[i]->pointLight.light_pos.getX(), spotlights[i]->pointLight.light_pos.getY(), spotlights[i]->pointLight.light_pos.getZ());
                //Point lightPosition = spotlights[i]->pointLight.light_pos;
                Point temp = (intersectionPoint - lightPosition).normalize();
                lightDirection.setValues(temp.getX(), temp.getY(), temp.getZ());

                double dotVal = lightDirection.dotProduct(spotlights[i]->light_direction);
                double tempVal = lightDirection.value() * spotlights[i]->light_direction.value();
                double angle = acos(dotVal / tempVal) * (180.0 / pi);

                if(fabs(angle) < spotlights[i]->cutoffAngle){

                    Ray lightRay(lightPosition, lightDirection);
                    Ray normal = getNormal(intersectionPoint, lightRay);
                    
                    double t2 = (intersectionPoint - lightPosition).value();
                    if(t2 < 1e-5) continue;
                    
                    //bool obscured = checkObscured(lightRay, t2, color);
                    
                    if(!checkObscured(lightRay, t2, color)){
                         handleObscured(color, colorAtIntersection, lightRay, normal , intersectionPoint , ray , i , false);  
                    }
                    else{
                        /// do nothing
                    
                    }
                }
            }
        }
        virtual int getNearestObjectIndex(Ray reflectionRay, vector<double> &color){
            int nearestObjectIndex = -1;
            double tMin = 1e9;
            for(int k = 0; k < objects.size(); k++){
                double t = objects[k]->intersect(reflectionRay, color, 0);
                if(t > 0 && t < tMin)
                    tMin = t , nearestObjectIndex = k;
            }

            return nearestObjectIndex;
        }   
        virtual void updateColorWithImpactOfReflection(vector<double> &color, vector<double> colorTemp , double reflectionCoefficient){
            color[0] += colorTemp[0] * reflectionCoefficient;
            color[1] += colorTemp[1] * reflectionCoefficient;
            color[2] += colorTemp[2] * reflectionCoefficient;
        }
		virtual double intersect(Ray ray, vector<double> &color, int level){
            double t = intersectHelper(ray, color, level);

            if(t < 0) return -1;
            if(level == 0) return t;
            Point tempPoint = ray.dir * t; 
            Point intersectionPoint = ray.start + tempPoint;
            auto colorAtIntersection = getColorAt(intersectionPoint);

            updateColorWithAmbience(color, colorAtIntersection , coefficients[0]);

            calculationsForPointLights(color, colorAtIntersection, intersectionPoint, ray);

            calculationsForSpotLights(color, colorAtIntersection, intersectionPoint, ray);

            if(level < recursionLevel){
             
                Ray normal = getNormal(intersectionPoint, ray);
                double tempVal = 2*(ray.dir.dotProduct(normal.dir));
                Point reflectionDir = ray.dir - normal.dir * tempVal;
                Ray reflectionRay(intersectionPoint, reflectionDir);
                Point normalizeReflecDir = reflectionDir.normalize();
                reflectionRay.start = reflectionRay.start + normalizeReflecDir * 1e-5;
                int nearestObjectIndex = getNearestObjectIndex(reflectionRay, color);

                if(nearestObjectIndex != -1){
                    vector<double> colorTemp(3,0); 
                    double t = objects[nearestObjectIndex]->intersect(reflectionRay,colorTemp, level+1);
                    updateColorWithImpactOfReflection(color, colorTemp , coefficients[3]);
                }
                
            }
            return t;
        }

        virtual ~Object(){
            coefficients.clear();
            coefficients.shrink_to_fit();
        }
};

class General : public Object{
    public:
    double A,B,C,D,E,F,G,H,I,J;

    General(double a , double b, double c, double d, double e, double f, double g, double h, double i, double j){
        A = a;
        B = b;
        C = c;
        D = d;
        E = e;
        F = f;
        G = g;
        H = h;
        I = i;
        J = j;
    }

    virtual void draw(){
        return;
    }

    virtual Ray getNormal(Point point, Ray incidentRay){
        Point temp;
        temp.x = 2 * A * point.x + D * point.y + E * point.z + G;
        temp.y = 2 * B * point.y + D * point.x + F * point.z + H;
        temp.z = 2 * C * point.z + E * point.x + F * point.y + I;
        Ray ray(point, temp);
        return ray;
    }

    bool isOk(Point point){
        return !((fabs(length)>1e-5 && (point.x < reference_point.x || point.x > reference_point.x + length)) || (fabs(width)>1e-5 && (point.y < reference_point.y || point.y > reference_point.y + width)) || (fabs(height)>1e-5 && (point.z < reference_point.z || point.z > reference_point.z + height)));
    }

    virtual void updateReturnStatus(double t1 , double t2 , Ray ray , pair<bool,double> &returnStatus){
        if(t1 > 0){
            Point intersectionPoint = ray.start + ray.dir * t1;
            if(isOk(intersectionPoint)){
                returnStatus = {true, t1};
            }
        }
        if(t2 > 0){
            Point intersectionPoint = ray.start + ray.dir * t2;
            if(isOk(intersectionPoint)){
                returnStatus = {true, t2};
            }
        }
    }

    virtual vector<double> calculate_c0_c1_c2(double x0, double y0, double z0, double x1, double y1, double z1){
        double c0 = A * pow(x1 , 2) + B * pow(y1 , 2) + C * pow(z1 , 2) + D * x1 * y1 + E * x1 * z1 + F * y1 * z1;
        double c1 = 2 * A * x0 * x1 + 2 * B * y0 * y1 + 2 * C * z0 * z1 + D * (x0 * y1 + x1 * y0) + E * (x0 * z1 + x1 * z0) + F * (y0 * z1 + y1 * z0) + G * x1 + H * y1 + I * z1;
        double c2 = A * pow(x0 , 2) + B * pow(y0 , 2) + C * pow(z0 , 2) + D * x0 * y0 + E * x0 * z0 + F * y0 * z0 + G * x0 + H * y0 + I * z0 + J;
        return {c0, c1, c2};
    }


    virtual double intersectHelper(Ray ray, vector<double> &color, int level){
        double x0 = ray.start.x;
        double y0 = ray.start.y;
        double z0 = ray.start.z;

        double x1 = ray.dir.x;
        double y1 = ray.dir.y;
        double z1 = ray.dir.z;
        double c0, c1, c2;

        auto list = calculate_c0_c1_c2(x0, y0, z0, x1, y1, z1);

        c0 = list[0];
        c1 = list[1];
        c2 = list[2];

        double discriminant = pow(c1 , 2) - 4 * c0 * c2;

        if(discriminant < 0) 
            return -1;

        if(fabs(c0) < 1e-5)
            return -(c2 / c1);

        double sqrtDiscriminant = sqrt(discriminant);
        
        double t1 = ((-1.0) * c1 - sqrtDiscriminant) / (2 * c0);
        double t2 = ((-1.0) * c1 + sqrtDiscriminant) / (2 * c0);

        if(t2 < t1) swap(t1,t2);

        pair<bool,double> returnStatus = {false, -1};

        updateReturnStatus(t1, t2, ray, returnStatus);

        return returnStatus.second;

    }

};


class Triangle: public Object{
    public:
    Point a , b , c;


    Triangle(Point a , Point b , Point c){
        this->a = a;
        this->b = b;
        this->c = c;
    }

    virtual Ray getReturnRay(Point normal , Ray incidentRay , Point point){
        Ray returnray(point, normal);
        if(incidentRay.dir.dotProduct(normal) < 0){
            returnray.dir = -normal;
        }
        
        return returnray;
    }

    virtual Ray getNormal(Point point, Ray incidentRay){
        Point normal = ((b-a).crossProduct(c-a)).normalize();

        Ray returnRay = getReturnRay(normal , incidentRay , point);

        return returnRay;
        
    }

    virtual void draw(){
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }

    virtual double determinant(double ar[3][3]){
	    double v1 = ar[0][0] * (ar[1][1] * ar[2][2] - ar[1][2] * ar[2][1]);
	    double v2 = (-1) * ar[0][1] * (ar[1][0] * ar[2][2] - ar[1][2] * ar[2][0]);
	    double v3 = ar[0][2] * (ar[1][0] * ar[2][1] - ar[1][1] * ar[2][0]);
        double det = v1 + v2 + v3;
	    return det;
    }

    virtual double intersectHelper(Ray ray, vector<double> &color, int level){
        double betaMatrix[3][3] = {{0}};
        double gammaMatrix[3][3] = {{0}};
        double tMatrix[3][3] = {{0}};
        double AMatrix[3][3] = {{0}};

        betaMatrix[0][0] = a.x - ray.start.x;
        betaMatrix[0][1] = a.x - c.x;  
        betaMatrix[0][2] = ray.dir.x;
        betaMatrix[1][0] = a.y - ray.start.y;
        betaMatrix[1][1] = a.y - c.y;
        betaMatrix[1][2] = ray.dir.y;
        betaMatrix[2][0] = a.z - ray.start.z;
        betaMatrix[2][1] = a.z - c.z;
        betaMatrix[2][2] = ray.dir.z;
        
        gammaMatrix[0][0] = a.x - b.x;
        gammaMatrix[0][1] = a.x - ray.start.x;
        gammaMatrix[0][2] = ray.dir.x;
        gammaMatrix[1][0] = a.y - b.y;
        gammaMatrix[1][1] = a.y - ray.start.y;
        gammaMatrix[1][2] = ray.dir.y;
        gammaMatrix[2][0] = a.z - b.z;
        gammaMatrix[2][1] = a.z - ray.start.z;
        gammaMatrix[2][2] = ray.dir.z;
		
        tMatrix[0][0] = a.x - b.x;
        tMatrix[0][1] = a.x - c.x;
        tMatrix[0][2] = a.x - ray.start.x;
        tMatrix[1][0] = a.y - b.y;
        tMatrix[1][1] = a.y - c.y;
        tMatrix[1][2] = a.y - ray.start.y;
        tMatrix[2][0] = a.z - b.z;
        tMatrix[2][1] = a.z - c.z;
        tMatrix[2][2] = a.z - ray.start.z;

        AMatrix[0][0] = a.x - b.x;
        AMatrix[0][1] = a.x - c.x;
        AMatrix[0][2] = ray.dir.x;
        AMatrix[1][0] = a.y - b.y;
        AMatrix[1][1] = a.y - c.y;
        AMatrix[1][2] = ray.dir.y;
        AMatrix[2][0] = a.z - b.z;
        AMatrix[2][1] = a.z - c.z;
        AMatrix[2][2] = ray.dir.z;

        double Adeterminant = determinant(AMatrix);
        double beta = determinant(betaMatrix) / Adeterminant;
        double gamma = determinant(gammaMatrix) / Adeterminant;
        double t = determinant(tMatrix) / Adeterminant;

        double returnVal = -1;

        if (beta > 0 && gamma > 0 && t > 0){
            if(beta + gamma < 1)
                returnVal = t;
            else{
                // do nothing
            }
        }

        return returnVal;
    }
};

class Sphere : public Object{
    public:

		Sphere(Point center, double radius){
			reference_point = center;
			length = radius;
		}

        virtual Ray getNormal(Point point, Ray incidentRay){
            Point temp = point - reference_point;
            return Ray(point, temp);
        }

        virtual void drawHemisphere(Point (&points)[100][100], int i, int j , int z){
            glBegin(GL_QUADS);{
                glVertex3f(points[i][j].x, points[i][j].y, z * points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, z * points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, z * points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, z * points[i + 1][j].z);
            }glEnd();
        }

        virtual void drawQuads(Point (&points)[100][100], int stacks, int slices){
            for (int i = 0; i < stacks; i++)
			{
				glPushMatrix();
				glTranslatef(this->reference_point.x, this->reference_point.y, this->reference_point.z);
				glColor3f(color[0], color[1], color[2]);
				for (int j = 0; j < slices; j++){
                    drawHemisphere(points , i , j , 1);
                    drawHemisphere(points , i , j , -1);
				}
				glPopMatrix();
			}    
        }

        virtual void generatePoints(Point (&points)[100][100], int stacks, int slices){
            for (int i = 0; i <= stacks; i++)
			{
				double h = length * sin(((double)i / (double)stacks) * (pi / 2));
				double r = length * cos(((double)i / (double)stacks) * (pi / 2));
				for (int j = 0; j <= slices; j++){
					points[i][j].x = r * cos(((double)j / (double)slices) * 2 * pi);
					points[i][j].y = r * sin(((double)j / (double)slices) * 2 * pi);
					points[i][j].z = h;
				}
			}
        }

		virtual void draw(){
            int stacks = 30;
			int slices = 20;

			Point points[100][100];
            generatePoints(points, stacks, slices);
            drawQuads(points, stacks, slices);
		}

        virtual double intersectHelper(Ray ray, vector<double> &color, int level){

            ray.start = ray.start - reference_point;
            
            double b = 2 * (ray.dir.dotProduct(ray.start));
            double c = (ray.start.dotProduct(ray.start)) - (length * length);

            double discriminant = pow(b, 2) - 4 * 1 * c;  /// a = 1
            double t = -1;

            double srqtDisc = sqrt(discriminant);
            
            double t1 = (-b - srqtDisc) / 2;
            double t2 = (-b + srqtDisc) / 2;

            if(t2 < t1) swap(t1, t2);

            if (t1 > 0 || t2 > 0){
                if(t1 > 0)
                    t = t1;
                else if (t2 > 0){
                    t = t2;
                }
            }
            
            return t;
        }

};

class Floor : public Object{
    public:
    int tiles;

    Floor(int floorWidth , int tileWidth){
        tiles = floorWidth / tileWidth;
        Point p(-floorWidth / 2, -floorWidth / 2, 0);
        reference_point = p;
        length = tileWidth;
    }

    virtual vector<double> getColorAt(Point point){

        int tileX = (point.x - reference_point.x) / length;
		int tileY = (point.y - reference_point.y) / length;
        if(tileX >= 0 && tileX < tiles && tileY >= 0 && tileY < tiles && (tileX + tileY) % 2 == 0){
            vector<double> color(3,1);
            return color;
        }
        else{
            vector<double> color(3,0);
            return color;
        }
    }

    virtual Ray getNormal(Point point, Ray incidentRay){
        Point temp(0, 0, 1);
        if(!(incidentRay.dir.z > 0)){
            temp.z = -1;
        }
        Ray ray(point, temp);
        return ray;
    }

    virtual void drawSquare(Point reference_point, double length, int i, int j){
        vector<double> ordinateValues(3);
        ordinateValues[0] = reference_point.x;
        ordinateValues[1] = reference_point.y;
        ordinateValues[2] = 0;
        glBegin(GL_QUADS);
		{
			glVertex3f(ordinateValues[0] + i * length, ordinateValues[1] + j * length, ordinateValues[2]);
			glVertex3f(ordinateValues[0] + (i + 1) * length, ordinateValues[1] + j * length, ordinateValues[2]);
			glVertex3f(ordinateValues[0] + (i + 1) * length, ordinateValues[1] + (j + 1) * length, ordinateValues[2]);
			glVertex3f(ordinateValues[0] + i * length, ordinateValues[1] + (j + 1) * length, ordinateValues[2]);
		}
		glEnd();
    }

    virtual void draw(){
        for (int i = 0; i < tiles; i++)
		{
			for (int j = 0; j < tiles; j++){
				if (((i + j) % 2) == 1){
                    glColor3f(0, 0, 0);
                }
				else{
                    glColor3f(1, 1, 1);
                }

                drawSquare(reference_point , length , i , j);
			}
		}
    }

    virtual double intersectHelper(Ray ray, vector<double> &color, int level){
        double rayDirz = ray.dir.z;
        double rayStartz = ray.start.z;
        
        if (round(rayDirz * 100) == 0)
			return -1;

        double t = -(rayStartz / rayDirz);

        return t;
    }
};