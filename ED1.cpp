/* 
 * File:   ED1.cpp
 * Author: Thomas Bronzwaer
 * Year:   2014
 * 
 * This C++ program creates a numerical plot of the Lienard-Wiechert potential
 * due to a moving charge (x(t) and v(t) specified by user).
 * 
 * The same field is also computed using a Fourier method (long-wavelength
 * approximation).
 * 
 * Units are cgs-Gauss. 
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

#include <cstdlib>
#include <complex>
#include "parameters_ED1.h"
#include "utilities.h"
#include "vector3.h"

using namespace std;

// DEFINITIONS
//////////////

#define _i_ complex<double>(0., 1.) // Imaginary unit
#define c   2.99792458e10           // Speed of light, cm/s

// FUNCTIONS
////////////

vector3 position(double t){
    vector3 X;
    X.x = AMPLITUDE * cos(FREQUENCY * t); // Horiz Osc
    //X.x = 0.5 * ACCELERATION * t * t; // Linear
    //X.y = SCALE * AMPLITUDE * sin(FREQUENCY * t); // Vert Osc
    //X.x = 0.;
    //X.y = VELOCITY * t;
    X.y = 0.;
    X.z = 0.;
    return X;
}

vector3 velocity(double t){
    vector3 V;
    V.x = -FREQUENCY * AMPLITUDE * sin(FREQUENCY * t); // Horiz Osc
    //V.x = ACCELERATION * t; // Linear
    //V.y = SCALE * FREQUENCY * AMPLITUDE * cos(FREQUENCY * t); // Vert Osc
    //V.x = 0.;
    //V.y = VELOCITY;
    V.y = 0.;
    V.z = 0.;
    return V;
}

// Returns the retarded time at the current position/time. 
// Uses binary search between 0 and t.
double ret_time(double t, vector3 r){
    double tau = t / 2.;   
    double divider = 4.;
    
    // Perform binary search NSTEPS times; accuracy ~t/2^NSTEPS
    for (int i = 2; i < NSTEPS; i++){
        if (tau < (t - separation(r, position(tau))/c)){
            tau += t / divider;
        }
        else{
            tau -= t / divider;
        }
        divider *= 2.;
    }
    
    return tau;
}

// Returns the advanced time at the current position/time. 
// Uses binary search between t and 2t ??????
double adv_time(double t, vector3 r){
    
    
    // DOESN'T WORK YET
    
    
    double tau = t * 1.5;   
    double divider = 4.;
    
    // Perform binary search NSTEPS times; accuracy ~t/2^NSTEPS
    for (int i = 2; i < NSTEPS; i++){
        if (tau < (t + separation(r, position(tau))/c)){
            tau += t / divider;
        }
        else{
            tau -= t / divider;
        }
        divider *= 2.;
    }
    
    return tau;
}

// MAIN FUNCTION
////////////////

int main() {
    
    // INITIALIZE VARIABLES
    ///////////////////////
    
    // RGB array for creating images
    unsigned char* data = new unsigned char[WIDTH * HEIGHT * 3];   
    FILE *file;
    file = fopen("C:/output/output.dat", "w");
    
    double k = FREQUENCY / c;           // Wave number
    vector3 xhat = vector3(1., 0., 0.);
    double X, Y, tau;
    vector3 r, R;
    double PHI, PHI_fourier;
    vector3 A, A_fourier, curl_A;
    int q = 0;                          // Image counter
    double RED,GRE,BLU;
    vector3 *A_field = new vector3[WIDTH * HEIGHT];
    
    // Set initial time
    double time = 0.0001;
    
    // MAIN COMPUTATION
    ///////////////////
    
    // for all time steps...
    for (int t = 0; t < TIMESTEPS; t++){
        
        // First pass over all pixels to compute the vector potential...
        for (int j = 0; j < HEIGHT; j++){
            for (int i = 0; i < WIDTH; i++){               
                // COMPUTE POSITION & RETARDED TIME
                ///////////////////////////////////

                // Position coordinates at each pixel center
                r.x = SCALE * ((double) i + .5 - WIDTH / 2.);
                //r.y = SCALE * ((double) j + .5 + 6. * HEIGHT / 1.);
                r.y = SCALE * ((double) j + .5 - HEIGHT / 2.);
                r.z = 0.;                
                
                tau = ret_time(time, r);               
                R = r - position(tau);      
                
                // COMPUTE LIÉNARD-WIECHERT SCALAR POTENTIAL (PHI)
                //////////////////////////////////////////////////

                PHI = CHARGE / (norm(R) - dot(R, velocity(tau))/c);

                // COMPUTE LIÉNARD-WIECHERT VECTOR POTENTIAL (A)
                ////////////////////////////////////////////////
                
                A = velocity(tau) / c * PHI;
                A_field[WIDTH * j + i] = A;
                
                // FOURIER VECTOR POTENTIAL (A_fourier)
                ///////////////////////////////////////
                
                complex<double> complexterm = (-_i_ * exp(_i_ * (k * norm(r) - 
                                              FREQUENCY * time)));
                double realpart = complexterm.real();
                A_fourier = xhat * k/norm(r) * CHARGE * AMPLITUDE * realpart;  
            }      
        }   
        
        // Second pass over all pixels - first we fill A, then we compute curl A
        for (int j = 0; j < HEIGHT; j++){
            for (int i = 0; i < WIDTH; i++){
                // COMPUTE B = curl A
                /////////////////////
                
                if (i > 1 && i < (WIDTH-1) && j > 1 && j < (HEIGHT-1)){
                    // Assume d/dz = 0 (we are in the XY-plane and the fields  
                    // are symmetric about this plane)
                    // We need the following derivatives;
                    double dAzdx = (A_field[WIDTH * j + (i+1)].z - 
                                    A_field[WIDTH * j + (i-1)].z) / (2.*SCALE);
                    double dAzdy = (A_field[WIDTH * (j+1) + i].z - 
                                    A_field[WIDTH * (j-1) + i].z) / (2.*SCALE);
                    double dAydx = (A_field[WIDTH * j + (i+1)].y - 
                                    A_field[WIDTH * j + (i-1)].y) / (2.*SCALE);
                    double dAxdy = (A_field[WIDTH * (j+1) + i].x - 
                                    A_field[WIDTH * (j-1) + i].x) / (2.*SCALE);
                    
                    curl_A.x = dAzdy - 0.;
                    curl_A.y = 0. - dAzdx;
                    curl_A.z = dAydx - dAxdy;
                
                }
                
                // COMPUTE PIXEL COLOR (PLOT)
                /////////////////////////////
                
                double factor = 12.e12;//6.e10;
                double color = norm(A);
                
                //RED = factor * color;
                //GRE = factor * color;
                //BLU = factor * color;
                
                //RED = factor * abs(A_field[WIDTH * j + i].x);
                //GRE = factor * abs(A_field[WIDTH * j + i].y);
                //BLU = factor * abs(A_field[WIDTH * j + i].z);
                
                //cout << "\n B.x " << (double) curl_A.x;
                //cout << "\n B.y " << (double) curl_A.y;
                //cout << "\n B.z " << (double) curl_A.z;
                
                RED = factor * abs(curl_A.x);
                GRE = factor * abs(curl_A.y);
                BLU = factor * abs(curl_A.z);

                RED = 255. * pow(RED, GAMMA);
                GRE = 255. * pow(GRE, GAMMA);
                BLU = 255. * pow(BLU, GAMMA);

                data[3 * (WIDTH * j + i) + 0] = min(RED,255.);
                data[3 * (WIDTH * j + i) + 1] = min(GRE,255.);
                data[3 * (WIDTH * j + i) + 2] = min(BLU,255.); 
                
                // PRINT DATA TO FILE
                /////////////////////
                
                if (j == 400){
                    fprintf(file, "\n%.8e %.8e", r.x, 
                            norm(A_field[WIDTH * j + i]));
                }
            }
        }
        
        write_image(data, q);
        q++;       
        time += T_INCREMENT;
    }

    // Cleanup & Exit
    /////////////////
    
    fclose(file);
    delete[] data;
    
    return 0;
}
