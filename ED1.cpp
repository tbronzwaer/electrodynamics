/* 
 * File:   ED1.cpp
 * Author: Thomas Bronzwaer
 * Year:   2014
 * 
 * This C++ program creates a numerical plot of the Lienard-Wiechert potential
 * due to a moving charge (x(t) and v(t) specified by user).
 * 
 * Units are cgs-Gauss. 
 */

#include <cstdlib>
#include "parameters_ED1.h"
#include "utilities.h"

using namespace std;

struct vector3{
    double x, y, z;
};

// FUNCTIONS
////////////

vector3 position(double t){
    vector3 X;
    X.x = SCALE * AMPLITUDE * cos(FREQUENCY * t); // Horiz Osc
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
    V.x = -SCALE * FREQUENCY * AMPLITUDE * sin(FREQUENCY * t); // Horiz Osc
    //V.x = ACCELERATION * t; // Linear
    //V.y = SCALE * FREQUENCY * AMPLITUDE * cos(FREQUENCY * t); // Vert Osc
    //V.x = 0.;
    //V.y = VELOCITY;
    V.y = 0.;
    V.z = 0.;
    return V;
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

// Returns the distance between the points A and B
double separation(vector3 A, vector3 B){
    vector3 C;
    C.x = A.x - B.x;
    C.y = A.y - B.y;
    C.z = A.z - B.z;
    return sqrt(dot(C, C));
}

vector3 A_to_B(vector3 A, vector3 B){   
    vector3 C;
    C.x = B.x - A.x;
    C.y = B.y - A.y;
    C.z = B.z - A.z;
    return C;
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

// MAIN FUNCTION
////////////////

int main() {
    
    // RGB array for creating images
    unsigned char* data = new unsigned char[WIDTH * HEIGHT * 3];
  
    // MAIN COMPUTATION
    ///////////////////
    
    double X, Y, tau;
    double RED,GRE,BLU;
    vector3 r, R;
    //double time = 1.5e-7;
    double time = 0.;
    double PHI = 0.;
    int q = 0; //img counter
    
    // for all times
    for (int t = 0; t < TIMESTEPS; t++){
        // For all pixels...
        for (int j = 0; j < HEIGHT; j++){
            for (int i = 0; i < WIDTH; i++){

                // Position coordinates at each pixel center
                r.x = SCALE * ((double) i + .5 - WIDTH / 2.);
                //r.y = SCALE * ((double) j + .5 + 6. * HEIGHT / 1.);
                r.y = SCALE * ((double) j + .5 - HEIGHT / 2.);
                r.z = 0.;                
                
                tau = ret_time(time, r);               
                R = A_to_B(position(tau), r);      
                   
                //PHI = 0.;               
                //if (sqrt(dot(R,R)) < c * time){
                if (1){
                    PHI = CHARGE / (sqrt(dot(R,R)) - dot(R, velocity(tau))/c);
                }
                
                // Set RGB at this location
                
                double factor = 7.;
                
                RED = factor * PHI;
                GRE = factor * PHI;
                BLU = factor * PHI;

                RED = 255. * pow(RED, GAMMA);
                GRE = 255. * pow(GRE, GAMMA);
                BLU = 255. * pow(BLU, GAMMA);

                data[3 * (WIDTH * j + i) + 0] = min(RED,255.);
                data[3 * (WIDTH * j + i) + 1] = min(GRE,255.);
                data[3 * (WIDTH * j + i) + 2] = min(BLU,255.); 
            }
        }
        
        write_image(data, q);
        q++;       
        time += T_INCREMENT;
    }

    // Cleanup & Exit
    /////////////////
    
    delete[] data;
    
    return 0;
}

