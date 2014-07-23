/* 
 * File:   parameters_ED1.h
 * Author: Thomas
 *
 * Created on July 7, 2014, 4:03 PM
 * 
 * NOTES ON cgs-Gauss UNITS:
 * 
 * Units are centimeter, gram, second
 */

#ifndef PARAMETERS_ED1_H
#define	PARAMETERS_ED1_H

#include <cstdlib>
#include <complex>
#include <stdio.h>
#include <iostream>

using namespace std;

// CONSTANTS
////////////

const double e               = 4.80320425e-10; // statcoulomb
const double c               = 2.99792458e10;  // cm/s

// PARAMETERS
/////////////

const int    WIDTH           = 800;
const int    HEIGHT          = 600;
const double SCALE           = 1.0; // Scale factor
const double VELOCITY        = 2.5e10;
const double ACCELERATION    = 4.e18;
const double FREQUENCY       = 3.e9;
const double AMPLITUDE       = 8.0;
const double CHARGE          = 1.0;
const int    NSTEPS          = 50; // For the binary search (ret. time)
const int    TIMESTEPS       = 1000;
const double T_INCREMENT     = 5.e-11;
const double GAMMA           = 1. / 2.2;
const std::string OUTPUT_DIR = "C:/output/electrodynamics";

#endif