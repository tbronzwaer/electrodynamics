/* 
 * File:   vector3.h (NOTE: only compatible with C++)
 * Author: Thomas Bronzwaer
 * 
 * This is a no-frills C++ implementation of ordinary Euclidean three-vectors, 
 * including the overloaded operators +, -, *, / to handle them more
 * conveniently as well as the norm, the dot- and cross products, and a function
 * that returns the distance between two points specified by vectors A & B.
 * 
 * Copyright (c) 2014 Thomas Bronzwaer
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy 
 * of this software and associated documentation files (the "Software"), to deal 
 * in the Software without restriction, including without limitation the rights 
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
 * copies of the Software, and to permit persons to whom the Software is 
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in 
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
 * SOFTWARE.
 */

#ifndef VECTOR3_H
#define	VECTOR3_H

#include <math.h>

class vector3{
    public:
        double x, y, z;
        vector3(){};
        vector3(double x_,double y_, double z_)
                :x(x_), y(y_), z(z_)
                {};
        double  dot       (const vector3& A, const vector3& B);
        double  norm      (const vector3& A);
        vector3 cross     (const vector3& A, const vector3& B);
        double  separation(const vector3& A, const vector3& B);
        vector3 operator +(const vector3 B); 
        vector3 operator -(const vector3 B); 
        vector3 operator *(const double a); 
        vector3 operator /(const double a); 
};

double dot(const vector3& A, const vector3& B){
    return A.x * B.x + A.y * B.y + A.z * B.z;
}

double norm(const vector3& A){
    return sqrt(dot(A,A));
}

vector3 cross(const vector3& A, const vector3& B){
    return vector3(A.y * B.z - A.z * B.y,
                   A.z * B.x - A.x * B.z,
                   A.x * B.y - A.y * B.x);
}

double separation(vector3 A, vector3 B){
    vector3 C;
    C.x = A.x - B.x;
    C.y = A.y - B.y;
    C.z = A.z - B.z;
    return norm(C);
}

vector3 vector3::operator+(const vector3 B){
    return vector3(x + B.x, y + B.y, z + B.z);
}

vector3 vector3::operator-(const vector3 B){
    return vector3(x - B.x, y - B.y, z - B.z);
}

vector3 vector3::operator*(const double a){
    return vector3(a*x, a*y, a*z);
}

vector3 vector3::operator/(const double a){
    return vector3(x/a, y/a, z/a);
}

#endif