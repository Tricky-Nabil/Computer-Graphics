#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<cmath>
#include<stack>
#include<algorithm>
#include<cstdlib>
#include <sys/stat.h>
#include "bitmap_image.hpp"

using namespace std;

const double PI = 2 * acos(0.0);

class Point{
public:
    double x , y , z , t;

    Point(double a , double b , double c) : x(a) , y(b) , z(c) , t(1){}

    Point() : t(1){}

    double dotProduct(Point p){
        double temp = x * p.x + y * p.y + z * p.z;
		return temp;
	}

    double value(){
		return sqrt(dotProduct(*this));
	}

	Point normalize(){
		return operator/(value());
	}

    Point crossProduct(Point p){
        Point temp;
        temp.x = y * p.z - z * p.y;
        temp.y = z * p.x - x * p.z;
        temp.z = x * p.y - y * p.x;

		return temp;
	}

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

    Point operator/(double val){
        Point temp;
        temp.x = x / val;
        temp.y = y / val;
        temp.z = z / val;

		return temp;
	}

    Point rotate(Point axis , double angle){
        Point dotPrdct = axis * this->dotProduct(axis);
        Point crossPrdct = axis.crossProduct(*this);
        double angleRad = angle * (PI / 180);

        return *this * cos(angleRad) + dotPrdct * (1 - cos(angleRad)) + crossPrdct * sin(angleRad);
    }
};

Point getUnitVector(int v){
        switch (v){
        case 1: return Point(1 , 0 , 0);
        case 2: return Point(0 , 1 , 0);
        case 3: return Point(0 , 0 , 1);
        default: return Point(0 , 0 , 0);
        }
        
}

class Triangle{
public:
    Point node1 , node2 , node3;

    Triangle(){}

    Triangle(Point n1, Point n2, Point n3) : node1(n1) , node2(n2) , node3(n3){}

};

double areaOfTr_2d(Triangle t){
        Point temp1 , temp2 , temp;
        temp1 = t.node2 - t.node1;
        temp2 = t.node3 - t.node1;
        temp = temp1.crossProduct(temp2);
        return temp.z;
}

Point getBaryCentricCo_ordinates(double x , double y , Triangle t){
        Point point(x , y , 0);
        double area = areaOfTr_2d(t);

        Triangle temp1(t.node2 , t.node3 , point) , temp2(t.node3 , t.node1 , point) , temp3(t.node1 , t.node2 , point);
        double a1 = areaOfTr_2d(temp1);
        double a2 = areaOfTr_2d(temp2);
        double a3 = areaOfTr_2d(temp3);

        Point temp;
        temp.x = a1 / area ;
        temp.y = a2 / area;
        temp.z = a3 / area;

        return temp;
}

vector<double> getBoundingBox(Triangle t){
        vector<double> v(4);
        double top , bottom , left , right;
        top = max({t.node1.y , t.node2.y , t.node3.y});
        bottom = min({t.node1.y , t.node2.y , t.node3.y});
        left = min({t.node1.x , t.node2.x , t.node3.x});
        right = max({t.node1.x , t.node2.x , t.node3.x});
        v[0] = top; v[1] = bottom; v[2] = left; v[3] = right;
        return v;

}

class Matrix{
    public:
    int row , column;
    vector<vector<double>> matrix;

    Matrix(int r = 4, int c = 4): row(r), column(c){
        matrix.resize(r , vector<double>(c));
    }

    void makeFirstColumn(Point p){
        matrix[0][0] = p.x;
        matrix[1][0] = p.y;
        matrix[2][0] = p.z;
        matrix[3][0] = 0;
    }

    void makeSecondColumn(Point p){
        matrix[0][1] = p.x;
        matrix[1][1] = p.y;
        matrix[2][1] = p.z;
        matrix[3][1] = 0;
    }

    void makeThirdColumn(Point p){
        matrix[0][2] = p.x;
        matrix[1][2] = p.y;
        matrix[2][2] = p.z;
        matrix[3][2] = 0;
    }

    void makeFourthColumn(){
        matrix[0][3] = 0;
        matrix[1][3] = 0;
        matrix[2][3] = 0;
        matrix[3][3] = 1;
    }

    void makeRotationMatrix(Point p1 , Point p2 , Point p3){
        makeFirstColumn(p1);
        makeSecondColumn(p2);
        makeThirdColumn(p3);
        makeFourthColumn();
    }

};

Point mul_Matrix_Point(Matrix m , Point point){
    Point temp(0 , 0 , 0);
    temp.t = 0;

    temp.x = m.matrix[0][0] * point.x + m.matrix[0][1] * point.y + m.matrix[0][2] * point.z + m.matrix[0][3] * point.t;
    temp.y = m.matrix[1][0] * point.x + m.matrix[1][1] * point.y + m.matrix[1][2] * point.z + m.matrix[1][3] * point.t;
    temp.z = m.matrix[2][0] * point.x + m.matrix[2][1] * point.y + m.matrix[2][2] * point.z + m.matrix[2][3] * point.t;
    temp.t = m.matrix[3][0] * point.x + m.matrix[3][1] * point.y + m.matrix[3][2] * point.z + m.matrix[3][3] * point.t;
        
    temp = temp / temp.t;
    return temp;
}

Triangle mul_Matrix_Triangle(Matrix m , Triangle t){
    Triangle temp;
    temp.node1 = mul_Matrix_Point(m , t.node1);
    temp.node2 = mul_Matrix_Point(m , t.node2);
    temp.node3 = mul_Matrix_Point(m , t.node3);
    return temp;
} 

Matrix mul_Matrix_Matrix(Matrix mat1 , Matrix mat2){
    Matrix mul;
    for(int i = 0; i < mat1.row ; i++){
        for(int j = 0 ; j < mat2.column ; j++){
            mul.matrix[i][j] = 0;
            for(int k = 0 ; k < mat1.column ; k++){
                mul.matrix[i][j] += mat1.matrix[i][k] * mat2.matrix[k][j];
            }
        }
    }
    return mul;
}

Matrix transpose(Matrix m){
        Matrix transposematrix;
        for(int i = 0 ; i < m.row ; i++){
            for(int j = 0; j < m.column ; j++){
                transposematrix.matrix[i][j] = m.matrix[j][i];
            }
        }
        return transposematrix;
}

Matrix getIdentityMatrix(){
        Matrix m;
        for(int i = 0 ; i < m.row ; i++){
            for (int j = 0 ; j < m.column ; j++){
                if(i == j) m.matrix[i][j] = 1;
                else m.matrix[i][j] = 0;
            }
        }
        return m;
}

Matrix getScaleMatrix(Point scale){
        Matrix m ;
        for(int i = 0 ; i < m.row ; i++){
            for (int j = 0 ; j < m.column ; j++){
                if(i == j) m.matrix[i][j] = 1;
                else m.matrix[i][j] = 0;
            }
        }
        m.matrix[0][0] = scale.x;
        m.matrix[1][1] = scale.y;
        m.matrix[2][2] = scale.z;
        return m;
}

Matrix getTranslateMatrix(Point translate){
        Matrix m ;
        for(int i = 0 ; i < m.row ; i++){
            for (int j = 0 ; j < m.column ; j++){
                if(i == j) m.matrix[i][j] = 1;
                else m.matrix[i][j] = 0;
            }
        }
        m.matrix[0][3] = translate.x;
        m.matrix[1][3] = translate.y;
        m.matrix[2][3] = translate.z;
        return m;
}

Matrix getRotationMatrix(Point rotate , double angle){
        rotate = rotate.normalize();

        Point c1 = getUnitVector(1).rotate(rotate , angle);
        Point c2 = getUnitVector(2).rotate(rotate , angle);
        Point c3 = getUnitVector(3).rotate(rotate , angle);

        Matrix rotationMatrix;
        rotationMatrix.makeRotationMatrix(c1 , c2 , c3);
        return rotationMatrix;
}

Matrix getProjectionMatrix(double fovY , double aspectRatio , double near , double far){
        double fovX = fovY * aspectRatio;
        double t = near * tan(fovY / 2);
        double r = near * tan(fovX / 2);

        Matrix m;
        for(int i = 0 ; i < m.row ; i++){
            for(int j = 0 ; j < m.column ; j++){
                m.matrix[i][j] = 0;
            }
        }
        m.matrix[0][0] = near / r;
        m.matrix[1][1] = near / t;
        m.matrix[2][2] = - (far + near) / (far - near);
        m.matrix[2][3] = - (2 * far * near) / (far - near);
        m.matrix[3][2] = -1;

        return m;
}

int main(int argc, char* argv[]){
    srand(1927);

    if(argc != 2){
        cout << "Two arguments required" << endl;
        return 1;
    }

    int result = mkdir("/home/nabil/OpenGL/Offline_2/output", 0777); // Create directory with permissions 777

    if (result == 0) {
       // std::cout << "Directory created successfully!" << std::endl;
    } else {
       // std::cerr << "Error creating directory: " << strerror(errno) << std::endl;
    }

    ifstream sceneInput(argv[1] + string("/scene.txt"));
    ifstream configInput(argv[1] + string("/config.txt"));

    ofstream stage1out("output/stage1.txt");
    ofstream stage2out("output/stage2.txt");
    ofstream stage3out("output/stage3.txt");
    ofstream zBufferoutput("output/z_buffer.txt");


    Point eye, look, up;
    double fovY, aspectRatio, near, far;
    double screenWidth, screenHeight;

    // reading from config file
    configInput >> screenWidth >> screenHeight;
    configInput.close();

    // reading from scene file
    sceneInput >>  eye.x >> eye.y >> eye.z;
    sceneInput >> look.x >> look.y >> look.z;
    sceneInput >> up.x >> up.y >> up.z;;
    sceneInput >> fovY >> aspectRatio >> near >> far;

    vector<Triangle> triangles;
    stack<Matrix> matrixStack;
    Matrix M = getIdentityMatrix();
    
    // stage1
    while(true) {
        string command;
        sceneInput >> command;

        if(command == "triangle"){
            Triangle t;
            sceneInput >> t.node1.x >> t.node1.y >> t.node1.z;
            sceneInput >> t.node2.x >> t.node2.y >> t.node2.z;
            sceneInput >> t.node3.x >> t.node3.y >> t.node3.z;

            t = mul_Matrix_Triangle(M , t);
            stage1out << fixed << showpoint << setprecision(7)  << t.node1.x << " " << setprecision(7) << t.node1.y << " " << setprecision(7) << t.node1.z<<" "<<endl;
            stage1out << fixed << showpoint << setprecision(7)  << t.node2.x << " " << setprecision(7) << t.node2.y << " " << setprecision(7) << t.node2.z<<" "<<endl;
            stage1out << fixed << showpoint << setprecision(7)  << t.node3.x << " " << setprecision(7) << t.node3.y << " " << setprecision(7) << t.node3.z<<" "<<endl;
            stage1out<<endl;
            triangles.push_back(t);
        }
        else if(command == "translate"){
            Point translateBy;
            sceneInput >> translateBy.x >> translateBy.y >> translateBy.z;

            M = mul_Matrix_Matrix(M , getTranslateMatrix(translateBy));
        }
        else if(command == "scale"){
            Point scaleBy;
            sceneInput >> scaleBy.x >> scaleBy.y >> scaleBy.z;

            M = mul_Matrix_Matrix(M , getScaleMatrix(scaleBy));
        }

        else if(command == "rotate"){
            double angle;
            Point rotateBy;
            sceneInput >> angle >> rotateBy.x >> rotateBy.y >> rotateBy.z;

            M = mul_Matrix_Matrix(M , getRotationMatrix(rotateBy , angle));
        }
        
        else if(command == "push"){
            matrixStack.push(M);
        }

        else if(command == "pop"){
            M = matrixStack.top();
            matrixStack.pop();
        }
        
        else if (command == "end"){
            break;
        }
    };
    
    stage1out.close();
    sceneInput.close();

  
    // stage 2
    Point l , r , u;
    l = (look - eye).normalize();
    r = l.crossProduct(up).normalize();
    u = r.crossProduct(l);

    Matrix eyeTranslate = getTranslateMatrix(-eye);
    Matrix temp;
    temp.makeRotationMatrix(r , u , l * (-1));

    Matrix viewRotate = transpose(temp); 
    
    Matrix viewTransformationMatrix = mul_Matrix_Matrix(viewRotate , eyeTranslate);

    for(auto &t : triangles){
        
        t = mul_Matrix_Triangle(viewTransformationMatrix , t);
        stage2out << fixed << showpoint << setprecision(7)  << t.node1.x << " " << setprecision(7) << t.node1.y << " " << setprecision(7) << t.node1.z<<" "<<endl;
        stage2out << fixed << showpoint << setprecision(7)  << t.node2.x << " " << setprecision(7) << t.node2.y << " " << setprecision(7) << t.node2.z<<" "<<endl;
        stage2out << fixed << showpoint << setprecision(7)  << t.node3.x << " " << setprecision(7) << t.node3.y << " " << setprecision(7) << t.node3.z<<" "<<endl;
        stage2out<<endl;
    }

    stage2out.close();

    
    // stage 3
    Matrix projectionTransformationMatrix = getProjectionMatrix(fovY * (PI / 180) , aspectRatio , near , far);

    // for(int i = 0 ; i < 4 ; i++){
    //     for (int j = 0 ; j < 4 ; j++){
    //         cout<<projectionTransformationMatrix.matrix[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    // int cnt = 0;
    for(auto &t : triangles){
        //t = projectionTransformationMatrix * t;
        // cout<<"triangle "<<cnt<<endl;
        // cout<<t.node1.x<<" "<<t.node1.y<<" "<<t.node1.z<<endl;
        // cout<<t.node2.x<<" "<<t.node2.y<<" "<<t.node2.z<<endl;
        // cout<<t.node3.x<<" "<<t.node3.y<<" "<<t.node3.z<<endl;
        // cout<<endl;
        // cnt++;
        t = mul_Matrix_Triangle(projectionTransformationMatrix , t);
        stage3out << fixed << showpoint << setprecision(7)  << t.node1.x << " " << setprecision(7) << t.node1.y << " " << setprecision(7) << t.node1.z<<" "<<endl;
        stage3out << fixed << showpoint << setprecision(7)  << t.node2.x << " " << setprecision(7) << t.node2.y << " " << setprecision(7) << t.node2.z<<" "<<endl;
        stage3out << fixed << showpoint << setprecision(7)  << t.node3.x << " " << setprecision(7) << t.node3.y << " " << setprecision(7) << t.node3.z<<" "<<endl;
        stage3out<<endl;
    }

    //stage 4
    vector<vector<double>> z_buffer;
    z_buffer.resize(screenHeight + 1 , vector<double>(screenWidth + 1 , 2));
    vector<vector<int>> frame_buffer;
    frame_buffer.resize(screenHeight + 1 , vector<int>(screenWidth + 1 , 49));

    double dx = 2 / screenWidth , dy = 2 / screenHeight;

    for(auto &t : triangles){
        int color = rand()% 30;

        double top , bottom , left , right;
        auto v = getBoundingBox(t);
        top = v[0] ; bottom = v[1] ; left = v[2] ; right = v[3];

        top = min(1.0 , top);
        bottom = max(-1.0 , bottom);
        left = max(-1.0 , left);
        right = min(1.0 , right);

        for(int y = top / dy  ; y >= bottom/dy - 1 ; y--){
            for(int x = (left / dx) - 1 ; x <= (right / dx) + 1 ; x++){
                Point temp = getBaryCentricCo_ordinates(x * dx , y * dy , t);
                double val = -1e-9;
                
                if(temp.x >= val && temp.y >= val && temp.z >= val){
                    Point temp = getBaryCentricCo_ordinates(x * dx , y * dy , t);
                    double z = temp.x * t.node1.z + temp.y * t.node2.z + temp.z * t.node3.z;

                    int row_index = screenHeight - (y + screenHeight / 2);
                    int col_index = x + screenWidth / 2;

                    if(row_index < 0 || row_index >= screenHeight || col_index < 0 || col_index >= screenWidth || z < -1 || z > 1) continue;
                    else if(z < z_buffer[row_index][col_index]){
                        z_buffer[row_index][col_index] = z;
                        frame_buffer[row_index][col_index] = color;
                    }
                }
            }
        }
    }

    //printing buffer

    for(int i = 0 ; i < z_buffer.size() ; i++){
        for(int j = 0 ; j < z_buffer[i].size() ; j++){
            if(z_buffer[i][j] < 1.1){
                zBufferoutput<<z_buffer[i][j]<<"\t";
            }
        }
        zBufferoutput<<endl;
    }

    zBufferoutput.close();

    // bitmap image generation
    bitmap_image image(screenWidth , screenHeight);
    for(int y = 0 ; y < screenHeight ; y++){
        for(int x = 0 ; x < screenWidth ; x++){
            int color = frame_buffer[y][x];
            image.set_pixel(x , y , palette_colormap[color]);
        }
    }
    image.save_image("output/out.bmp");
}