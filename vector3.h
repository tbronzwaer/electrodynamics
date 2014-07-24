/* 
 * File:   vector3.h
 * Author: Thomas
 *
 * Created on July 24, 2014, 11:01 AM
 */

#ifndef VECTOR3_H
#define	VECTOR3_H

class vector3{
    public:
        double x, y, z;
    vector3 operator+(const vector3 B); 
    vector3 operator-(const vector3 B); 
    vector3 operator*(const double a); 
    vector3 operator/(const double a); 
};

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

double dot(const vector3& A, const vector3& B){
    return A.x * B.x + A.y * B.y + A.z * B.z;
}

vector3 cross(const vector3& A, const vector3& B){
    vector3 C;
    C.x = A.y * B.z - A.z * B.y;
    C.y = A.z * B.x - A.x * B.z;
    C.z = A.x * B.y - A.y * B.x;
    return C;
}

#endif	/* VECTOR3_H */

