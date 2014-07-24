/* 
 * File:   vector3.h
 * Author: Thomas Bronzwaer
 *
 * Created on July 24, 2014, 11:01 AM
 * 
 * This is a fairly 'vanilla' implementation of ordinary Euclidean three-
 * vectors, including the overloaded operators +, -, *, / to handle them more
 * conveniently as well as the norm, the dot- and cross products, and a function
 * that returns the distance between two points specified by vectors A & B.
 */

#ifndef VECTOR3_H
#define	VECTOR3_H

class vector3{
    public:
        double x, y, z;
        vector3(){};
        vector3(double x_,double y_, double z_)
                :x(x_), y(y_), z(z_)
                {};
        double dot(const vector3& A, const vector3& B);
        double norm(const vector3& A);
        vector3 cross(const vector3& A, const vector3& B);
        double separation(const vector3& A, const vector3& B);
        vector3 operator+(const vector3 B); 
        vector3 operator-(const vector3 B); 
        vector3 operator*(const double a); 
        vector3 operator/(const double a); 
};


double dot(const vector3& A, const vector3& B){
    return A.x * B.x + A.y * B.y + A.z * B.z;
}

double norm(const vector3& A){
    return sqrt(dot(A,A));
}

vector3 cross(const vector3& A, const vector3& B){
    vector3 C;
    C.x = A.y * B.z - A.z * B.y;
    C.y = A.z * B.x - A.x * B.z;
    C.z = A.x * B.y - A.y * B.x;
    return C;
}

// Returns the distance between the points A and B
double separation(vector3 A, vector3 B){
    vector3 C;
    C.x = A.x - B.x;
    C.y = A.y - B.y;
    C.z = A.z - B.z;
    return sqrt(dot(C, C));
}

vector3 vector3::operator+(const vector3 B){
    vector3 C;
    C.x = x + B.x;
    C.y = y + B.y;
    C.z = z + B.z;
    return C;
}

vector3 vector3::operator-(const vector3 B){
    vector3 C;
    C.x = x - B.x;
    C.y = y - B.y;
    C.z = z - B.z;
    return C;
}

vector3 vector3::operator*(const double a){
    vector3 C;
    C.x = x * a;
    C.y = y * a;
    C.z = z * a;
    return C;
}

vector3 vector3::operator/(const double a){
    vector3 C;
    C.x = x / a;
    C.y = y / a;
    C.z = z / a;
    return C;
}

#endif