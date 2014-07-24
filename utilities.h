/* 
 * File:   utilities.h
 * Author: Thomas
 *
 * Created on March 21, 2014, 12:43 AM
 * 
 * This header file contains the following utility functions for 3DPauli:
 * 
 * Complex to RGB conversion
 * BMP image writing routine for Windows
 */

#ifndef UTILITIES_H
#define	UTILITIES_H

#include <windows.h>
#include <sstream>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "parameters_ED1.h"


using namespace std;

void complex_to_rgb(const complex<double> C, double RGB[]){
    // PLOT THE COMPLEX NUMBER C  
    double magnitude = abs(C);
    double phase = arg(C);

    // Normalize the phase.
    phase += M_PI;
    phase /= (2. * M_PI);

    // Normalize the magnitude.
    magnitude /= 1.0;

    // HSV to RGB method from Wikipedia 
    double Hprime = phase * 6.;	
    double Ch = magnitude * 1.;	
    double X = Ch * (1. - abs(fmod(Hprime, 2.) - 1.));

    X = max(0., min(255., X));
    Ch = max(0., min(255., Ch));

    if (0.0 <= Hprime && Hprime <= 1.0)
        { RGB[0] = Ch; RGB[1] = X; RGB[2] = 0.; }
    else if (1.0 <= Hprime && Hprime <= 2.0)
        { RGB[0] = X; RGB[1] = Ch; RGB[2] = 0.; }
    else if (2.0 <= Hprime && Hprime <= 3.0)
        { RGB[0] = 0.; RGB[1] = Ch; RGB[2] = X; }
    else if (3.0 <= Hprime && Hprime <= 4.0)
        { RGB[0] = 0.; RGB[1] = X; RGB[2] = Ch; }
    else if (4.0 <= Hprime && Hprime <= 5.0)
        { RGB[0] = X; RGB[1] = 0.; RGB[2] = Ch; }
    else if (5.0 <= Hprime && Hprime <= 6.0)
        { RGB[0] = Ch; RGB[1] = 0.; RGB[2] = X; }
    else 
        { RGB[0] = 0.; RGB[1] = 0.; RGB[2] = 0.; }
    // END <- material from Wikipedia
}

int write_image(unsigned char *data, int imagecounter){
    stringstream filename;
    string filenamestring;
    filename.str("");
    if (imagecounter<10)
        filename << OUTPUT_DIR << "/frame0000" << (int) imagecounter << ".bmp";
    else if (imagecounter<100)
        filename << OUTPUT_DIR << "/frame000" << (int) imagecounter << ".bmp";
    else if (imagecounter<1000)
        filename <<  OUTPUT_DIR << "/frame00" << (int) imagecounter << ".bmp";
    else if (imagecounter<10000)
        filename <<  OUTPUT_DIR << "/frame0" << (int) imagecounter << ".bmp";
    else
        filename << OUTPUT_DIR << "/frame" << (int) imagecounter << ".bmp";

    filenamestring = filename.str();

    HANDLE file;
    BITMAPFILEHEADER header;
    BITMAPINFOHEADER info;
    RGBTRIPLE *image;
    DWORD write = 0;

    image = new RGBTRIPLE[WIDTH * HEIGHT];
    file = CreateFile(filenamestring.c_str(), GENERIC_WRITE, 0, NULL, 
                      CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL); 
    //ofstream file2;
    //file2.open(filenamestring.c_str());
    //file2.close();

    header.bfType = 19778;                                                                    //Sets our type to BM or bmp
    header.bfSize = sizeof(header.bfOffBits) + sizeof(RGBTRIPLE);                                                //Sets the size equal to the size of the header struct
    header.bfReserved1 = 0;                                                                    //sets the reserves to 0
    header.bfReserved2 = 0;
    header.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);                    //Sets offbits equal to the size of file and info header

    info.biSize = sizeof(BITMAPINFOHEADER);
    info.biWidth = WIDTH;
    info.biHeight = HEIGHT;
    info.biPlanes = 1;
    info.biBitCount = 24;
    info.biCompression = BI_RGB;
    info.biSizeImage = WIDTH * HEIGHT * (24/8);
    info.biXPelsPerMeter = 2400;
    info.biYPelsPerMeter = 2400;
    info.biClrImportant = 0;
    info.biClrUsed = 0;

    WriteFile(file, &header, sizeof(header), &write, NULL);
    WriteFile(file, &info, sizeof(info), &write, NULL);

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            image[WIDTH * y + x].rgbtRed = data[(WIDTH*y+x)*3 + 0];
            image[WIDTH * y + x].rgbtGreen = data[(WIDTH*y+x)*3+1];
            image[WIDTH * y + x].rgbtBlue = data[(WIDTH*y+x)*3 +2];
        }
    }
    WriteFile(file, image, info.biSizeImage, &write, NULL);
    CloseHandle(file);    

    return 0;
}

#endif