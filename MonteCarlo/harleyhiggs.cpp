// Created 12 July 2019
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <complex.h>

using namespace std;
using namespace std::chrono;


#define I1 complex <double> (0., 1.) // imaginary unit
#define PI 3.141592653589793238462643383279502884197169399375105820974944



inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}


long double gXi, PXi, b1, b2, ac;                           // NLsigModel coupling
int Lx, Ly, Lt;                         // number of spins in x and y and t (imaginary time)
int N;                                  // number of spins




double T, T0, H, Ps;                               // temperature

long double theta0, thetaD, alpha_c;                  // angles
long double phi0, phiD;






//Rotations
double  phiH, thetaH, AmpH, phiH1, thetaH1, AmpH1, phiH2, thetaH2, AmpH2, phiH3, thetaH3, AmpH3, phiH4, thetaH4, AmpH4;
double d_Phi, alpha, C_alpha, S_alpha;
double d_Amp0, d_Amp, d_A;
double thetaAx, thetaAy, phiAx, phiAy;

//double  phiAx, thetaAx, phiAy, thetaAy;



//Coupling Constants
double ps, ms, u0, u1, u2, u3, u4, g, Keff, Geff;

//Higgs and Gauge Amplitude
double H0;
double A0x, A0y;
int a; //random integer: 0 or 1
double rH;



//Defect Locations
int s1, s2, s3, s4, s5, k1, k2, k3, k4, k5;
int Adefect, Hdefect, Ldefect;


//----- Define variables -----

double acceptanceRatioH1;
double acceptanceRatioH1Rab;
double acceptanceRatioH1Rbc;
double acceptanceRatioH1Rac;
double acceptanceRatioH2;
double acceptanceRatioH2Rab;
double acceptanceRatioH2Rbc;
double acceptanceRatioH2Rac;
double acceptanceRatioH3;
double acceptanceRatioH3Rab;
double acceptanceRatioH3Rbc;
double acceptanceRatioH3Rac;
double acceptanceRatioH4;
double acceptanceRatioH4Rab;
double acceptanceRatioH4Rbc;
double acceptanceRatioH4Rac;
double XIRatio;
double acceptanceRatioAx;
double acceptanceRatioAxRab;
double acceptanceRatioAxRbc;
double acceptanceRatioAxRac;
double acceptanceRatioAy;
double acceptanceRatioAyRab;
double acceptanceRatioAyRbc;
double acceptanceRatioAyRac;

double acceptsAx = 0.;
double acceptsAy = 0.;
double acceptsAxAmp = 0.;
double acceptsAyAmp = 0.;
//--------------------------------


//Fixed global angles

double comrandH1x;
double comrandH1y;
double comrandH2x;
double comrandH2y;
double comrandH3x;
double comrandH3y;
double comrandH4x;
double comrandH4y;


//            comrandH1x = std_rand();
//            comrandH1y = std_rand();
//            comrandH2x = std_rand();
//            comrandH2y = std_rand();
//            comrandH3x = std_rand();
//            comrandH3y = std_rand();
//            comrandH4x = std_rand();
//            comrandH4y = std_rand();



int steps = 0;                          // steps so far



using std::vector;


vector< vector< vector<double> > >   H1a;
vector< vector< vector<double> > >   H1b;
vector< vector< vector<double> > >   H1c;


vector< vector< vector<double> > >  H2a;
vector< vector< vector<double> > >  H2b;
vector< vector< vector<double> > >  H2c;


vector< vector< vector<double> > >  H3a;
vector< vector< vector<double> > >  H3b;
vector< vector< vector<double> > >  H3c;


vector< vector< vector<double> > >  H4a;
vector< vector< vector<double> > >  H4b;
vector< vector< vector<double> > >  H4c;





vector< vector< vector<double> > >  H1a0;
vector< vector< vector<double> > >  H1b0;
vector< vector< vector<double> > >  H1c0;

vector< vector< vector<double> > >  H2a0;
vector< vector< vector<double> > >  H2b0;
vector< vector< vector<double> > >  H2c0;

vector< vector< vector<double> > >  H3a0;
vector< vector< vector<double> > >  H3b0;
vector< vector< vector<double> > >  H3c0;

vector< vector< vector<double> > >  H4a0;
vector< vector< vector<double> > >  H4b0;
vector< vector< vector<double> > >  H4c0;





vector< vector< vector<double> > >  H1aflip;
vector< vector< vector<double> > >  H1bflip;
vector< vector< vector<double> > >  H1cflip;

vector< vector< vector<double> > >  H2aflip;
vector< vector< vector<double> > >  H2bflip;
vector< vector< vector<double> > >  H2cflip;

vector< vector< vector<double> > >  H3aflip;
vector< vector< vector<double> > >  H3bflip;
vector< vector< vector<double> > >  H3cflip;

vector< vector< vector<double> > >  H4aflip;
vector< vector< vector<double> > >  H4bflip;
vector< vector< vector<double> > >  H4cflip;






void InitializeConfiguration () {
    
    
    H1a.resize(Lx);
    H1b.resize(Lx);
    H1c.resize(Lx);
    
    H2a.resize(Lx);
    H2b.resize(Lx);
    H2c.resize(Lx);
    
    H3a.resize(Lx);
    H3b.resize(Lx);
    H3c.resize(Lx);
    
    H4a.resize(Lx);
    H4b.resize(Lx);
    H4c.resize(Lx);

    
    
    
    H1a0.resize(Lx);
    H1b0.resize(Lx);
    H1c0.resize(Lx);
    
    H2a0.resize(Lx);
    H2b0.resize(Lx);
    H2c0.resize(Lx);
    
    H3a0.resize(Lx);
    H3b0.resize(Lx);
    H3c0.resize(Lx);
    
    H4a0.resize(Lx);
    H4b0.resize(Lx);
    H4c0.resize(Lx);
    

    
    
    
    H1aflip.resize(Lx);
    H1bflip.resize(Lx);
    H1cflip.resize(Lx);
    
    H2aflip.resize(Lx);
    H2bflip.resize(Lx);
    H2cflip.resize(Lx);
    
    H3aflip.resize(Lx);
    H3bflip.resize(Lx);
    H3cflip.resize(Lx);
    
    H4aflip.resize(Lx);
    H4bflip.resize(Lx);
    H4cflip.resize(Lx);
    

    
    
    for (int i = 0; i < Lx; i++){


        
        H1a[i].resize(Ly);
        H1b[i].resize(Ly);
        H1c[i].resize(Ly);
        
        H2a[i].resize(Ly);
        H2b[i].resize(Ly);
        H2c[i].resize(Ly);
        
        H3a[i].resize(Ly);
        H3b[i].resize(Ly);
        H3c[i].resize(Ly);
        
        H4a[i].resize(Ly);
        H4b[i].resize(Ly);
        H4c[i].resize(Ly);
        
        H1a0[i].resize(Ly);
        H1b0[i].resize(Ly);
        H1c0[i].resize(Ly);
        
        H2a0[i].resize(Ly);
        H2b0[i].resize(Ly);
        H2c0[i].resize(Ly);
        
        H3a0[i].resize(Ly);
        H3b0[i].resize(Ly);
        H3c0[i].resize(Ly);
        
        H4a0[i].resize(Ly);
        H4b0[i].resize(Ly);
        H4c0[i].resize(Ly);
        
        H1aflip[i].resize(Ly);
        H1bflip[i].resize(Ly);
        H1cflip[i].resize(Ly);
        
        H2aflip[i].resize(Ly);
        H2bflip[i].resize(Ly);
        H2cflip[i].resize(Ly);
        
        H3aflip[i].resize(Ly);
        H3bflip[i].resize(Ly);
        H3cflip[i].resize(Ly);
        
        H4aflip[i].resize(Ly);
        H4bflip[i].resize(Ly);
        H4cflip[i].resize(Ly);
        

        
        
    }
    
    
    
    for (int i = 0; i < Lx; i++){
        for (int j = 0; j < Ly; j++){
        
        
        
        H1a[i][j].resize(Lt);
        H1b[i][j].resize(Lt);
        H1c[i][j].resize(Lt);
        
        H2a[i][j].resize(Lt);
        H2b[i][j].resize(Lt);
        H2c[i][j].resize(Lt);
        
        H3a[i][j].resize(Lt);
        H3b[i][j].resize(Lt);
        H3c[i][j].resize(Lt);
        
        H4a[i][j].resize(Lt);
        H4b[i][j].resize(Lt);
        H4c[i][j].resize(Lt);
        
        H1a0[i][j].resize(Lt);
        H1b0[i][j].resize(Lt);
        H1c0[i][j].resize(Lt);
        
        H2a0[i][j].resize(Lt);
        H2b0[i][j].resize(Lt);
        H2c0[i][j].resize(Lt);
        
        H3a0[i][j].resize(Lt);
        H3b0[i][j].resize(Lt);
        H3c0[i][j].resize(Lt);
        
        H4a0[i][j].resize(Lt);
        H4b0[i][j].resize(Lt);
        H4c0[i][j].resize(Lt);
        
        H1aflip[i][j].resize(Lt);
        H1bflip[i][j].resize(Lt);
        H1cflip[i][j].resize(Lt);
        
        H2aflip[i][j].resize(Lt);
        H2bflip[i][j].resize(Lt);
        H2cflip[i][j].resize(Lt);
        
        H3aflip[i][j].resize(Lt);
        H3bflip[i][j].resize(Lt);
        H3cflip[i][j].resize(Lt);
        
        H4aflip[i][j].resize(Lt);
        H4bflip[i][j].resize(Lt);
        H4cflip[i][j].resize(Lt);
        
        
        
        
    }
    }
    

    
    
    for (int i = 0; i < Lx; i++){
        for (int j = 0; j < Ly; j++){
            for (int t = 0; t < Lt; t++){
            

            
            // Initialize Higgs Fields
            
            //Local angles
            comrandH1x = std_rand();
            comrandH1y = std_rand();
            comrandH2x = std_rand();
            comrandH2y = std_rand();
            comrandH3x = std_rand();
            comrandH3y = std_rand();
            comrandH4x = std_rand();
            comrandH4y = std_rand();

            
//            double rH = 1.0; //WARNING 1.0 for F-G Phase and 0.1 for C Phase; //rand()%2;
            
              phiH =2.*PI*comrandH1x; //0.*(2.-2.*std_rand())*PI;                             // generate random angles for Higgs H1
              thetaH = PI*comrandH1y; //0.*(2.-2.*std_rand())*PI;                          // generate random angles for Higgs H1
              AmpH = H0*1.+d_Amp0*(1.-2.*std_rand());   //WARNING: No rH*             // generate random amplitudes for Higgs H1
            
              phiH2 =2.*PI*comrandH2x; //0.*(2.-2.*std_rand())*PI;                             // generate random angles for Higgs H1
              thetaH2 = PI*comrandH2y; //0.*(2.-2.*std_rand())*PI;                          // generate random angles for Higgs H1
              AmpH2 = H0*1.+d_Amp0*(1.-2.*std_rand());   //WARNING: No rH*           // generate random amplitudes for Higgs H1
            
              phiH3 =2.*PI*comrandH3x;                             // generate random angles for Higgs H1
              thetaH3 = PI*comrandH3y;                          // generate random angles for Higgs H1
              AmpH3 = H0*1.+d_Amp0*(1.-2.*std_rand());             // generate random amplitudes for Higgs H1
            
              phiH4 =2.*PI*comrandH4x;                             // generate random angles for Higgs H1
              thetaH4 = PI*comrandH4y;                          // generate random angles for Higgs H1
              AmpH4 = H0*1.+d_Amp0*(1.-2.*std_rand());             // generate random amplitudes for Higgs H1
            
            
            H1a[i][j][t] =double(AmpH*sin(thetaH)*cos(phiH));
            H1b[i][j][t] =double(AmpH*sin(thetaH)*sin(phiH));
            H1c[i][j][t] =double(AmpH*cos(thetaH));
            
            H2a[i][j][t] =double(AmpH2*sin(thetaH2)*cos(phiH2));
            H2b[i][j][t] =double(AmpH2*sin(thetaH2)*sin(phiH2));
            H2c[i][j][t] =double(AmpH2*cos(thetaH2));
            
            H3a[i][j][t] =double(AmpH3*sin(thetaH3)*cos(phiH3));
            H3b[i][j][t] =double(AmpH3*sin(thetaH3)*sin(phiH3));
            H3c[i][j][t] =double(AmpH3*cos(thetaH3));
            
            H4a[i][j][t] =double(AmpH4*sin(thetaH4)*cos(phiH4));
            H4b[i][j][t] =double(AmpH4*sin(thetaH4)*sin(phiH4));
            H4c[i][j][t] =double(AmpH4*cos(thetaH4));
            
            

            }}}
        
        
}


std::complex<double> EConfiguration(int s, int k, int w){
    

    std::complex<double> VH=0.+0.*I1;
    std::complex<double> DiH=0.+0.*I1;
    std::complex<double> DjH=0.+0.*I1;
    std::complex<double> DtH=0.+0.*I1;
    
    std::complex<double> DijHH=0.+0.*I1;
    //std::complex<double> DitHH=0.+0.*I1;
    //std::complex<double> DjiHH=0.+0.*I1;
    //std::complex<double> DtiHH=0.+0.*I1;
    
    //std::complex<double> DjtHH=0.+0.*I1;
    std::complex<double> DtjHH=0.+0.*I1;

    
    std::complex<double> eConfig=0.+0.*I1;
   
    
    
    int sPrev = s == 0 ? Lx-1 : s-1;
    int sNext = s == Lx-1 ? 0 : s+1;
    int kPrev = k == 0 ? Ly-1 : k-1;
    int kNext = k == Ly-1 ? 0 : k+1;
    int wPrev = w == 0 ? Lt-1 : w-1;
    int wNext = w == Lt-1 ? 0 : w+1;
    
    
    //Higgs quadratic and Quartic Potential
    
    VH = ms*(pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2) +
             pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2) +
             pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2) +
             pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
    u1*(((H1a[s][k][w] + I1*H2a[s][k][w])*
         (H3a[s][k][w] - I1*H4a[s][k][w]) +
         (H1b[s][k][w] + I1*H2b[s][k][w])*
         (H3b[s][k][w] - I1*H4b[s][k][w]) +
         (H1c[s][k][w] + I1*H2c[s][k][w])*
         (H3c[s][k][w] - I1*H4c[s][k][w]))*
        ((H1a[s][k][w] - I1*H2a[s][k][w])*
         (H3a[s][k][w] + I1*H4a[s][k][w]) +
         (H1b[s][k][w] - I1*H2b[s][k][w])*
         (H3b[s][k][w] + I1*H4b[s][k][w]) +
         (H1c[s][k][w] - I1*H2c[s][k][w])*
         (H3c[s][k][w] + I1*H4c[s][k][w])) +
        ((H1a[s][k][w] - I1*H2a[s][k][w])*
         (H3a[s][k][w] - I1*H4a[s][k][w]) +
         (H1b[s][k][w] - I1*H2b[s][k][w])*
         (H3b[s][k][w] - I1*H4b[s][k][w]) +
         (H1c[s][k][w] - I1*H2c[s][k][w])*
         (H3c[s][k][w] - I1*H4c[s][k][w]))*
        ((H1a[s][k][w] + I1*H2a[s][k][w])*
         (H3a[s][k][w] + I1*H4a[s][k][w]) +
         (H1b[s][k][w] + I1*H2b[s][k][w])*
         (H3b[s][k][w] + I1*H4b[s][k][w]) +
         (H1c[s][k][w] + I1*H2c[s][k][w])*
         (H3c[s][k][w] + I1*H4c[s][k][w])) +
        ((pow(H1a[s][k][w] - I1*H2a[s][k][w],2) +
          pow(H1b[s][k][w] - I1*H2b[s][k][w],2) +
          pow(H1c[s][k][w] - I1*H2c[s][k][w],2))*
         (pow(H1a[s][k][w] + I1*H2a[s][k][w],2) +
          pow(H1b[s][k][w] + I1*H2b[s][k][w],2) +
          pow(H1c[s][k][w] + I1*H2c[s][k][w],2)) +
         (pow(H3a[s][k][w] - I1*H4a[s][k][w],2) +
          pow(H3b[s][k][w] - I1*H4b[s][k][w],2) +
          pow(H3c[s][k][w] - I1*H4c[s][k][w],2))*
         (pow(H3a[s][k][w] + I1*H4a[s][k][w],2) +
          pow(H3b[s][k][w] + I1*H4b[s][k][w],2) +
          pow(H3c[s][k][w] + I1*H4c[s][k][w],2)))/2. +
        pow(pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2) +
            pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2) -
            pow(H3a[s][k][w],2) - pow(H3b[s][k][w],2) - pow(H3c[s][k][w],2) -
            pow(H4a[s][k][w],2) - pow(H4b[s][k][w],2) - pow(H4c[s][k][w],2),2)/4.) +
    u0*pow(pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2) +
           pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2) +
           pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2) +
           pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2),2);
    


    DiH = Keff*(2*(2*(H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] +
                      H1c[s][k][w]*H2c[s][k][w])*
                   (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] +
                    H1c[sNext][k][w]*H2c[sNext][k][w]) +
                   2*(H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] +
                      H1c[s][k][w]*H3c[s][k][w])*
                   (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] +
                    H1c[sNext][k][w]*H3c[sNext][k][w]) +
                   2*(H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] +
                      H2c[s][k][w]*H3c[s][k][w])*
                   (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] +
                    H2c[sNext][k][w]*H3c[sNext][k][w]) +
                   2*(H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] +
                      H1c[s][k][w]*H4c[s][k][w])*
                   (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] +
                    H1c[sNext][k][w]*H4c[sNext][k][w]) +
                   2*(H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] +
                      H2c[s][k][w]*H4c[s][k][w])*
                   (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] +
                    H2c[sNext][k][w]*H4c[sNext][k][w]) +
                   2*(H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] +
                      H3c[s][k][w]*H4c[s][k][w])*
                   (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] +
                    H3c[sNext][k][w]*H4c[sNext][k][w]) +
                   (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*
                   (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2))\
                   + (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*
                   (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2))\
                   + (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*
                   (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2))\
                   + (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*
                   (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2))))/3.;
    
    DjH =  Keff*(2*(2*(H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] +
                       H1c[s][k][w]*H2c[s][k][w])*
                    (H1a[s][kNext][w]*H2a[s][kNext][w] + H1b[s][kNext][w]*H2b[s][kNext][w] +
                     H1c[s][kNext][w]*H2c[s][kNext][w]) +
                    2*(H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] +
                       H1c[s][k][w]*H3c[s][k][w])*
                    (H1a[s][kNext][w]*H3a[s][kNext][w] + H1b[s][kNext][w]*H3b[s][kNext][w] +
                     H1c[s][kNext][w]*H3c[s][kNext][w]) +
                    2*(H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] +
                       H2c[s][k][w]*H3c[s][k][w])*
                    (H2a[s][kNext][w]*H3a[s][kNext][w] + H2b[s][kNext][w]*H3b[s][kNext][w] +
                     H2c[s][kNext][w]*H3c[s][kNext][w]) +
                    2*(H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] +
                       H1c[s][k][w]*H4c[s][k][w])*
                    (H1a[s][kNext][w]*H4a[s][kNext][w] + H1b[s][kNext][w]*H4b[s][kNext][w] +
                     H1c[s][kNext][w]*H4c[s][kNext][w]) +
                    2*(H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] +
                       H2c[s][k][w]*H4c[s][k][w])*
                    (H2a[s][kNext][w]*H4a[s][kNext][w] + H2b[s][kNext][w]*H4b[s][kNext][w] +
                     H2c[s][kNext][w]*H4c[s][kNext][w]) +
                    2*(H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] +
                       H3c[s][k][w]*H4c[s][k][w])*
                    (H3a[s][kNext][w]*H4a[s][kNext][w] + H3b[s][kNext][w]*H4b[s][kNext][w] +
                     H3c[s][kNext][w]*H4c[s][kNext][w]) +
                    (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*
                    (pow(H1a[s][kNext][w],2) + pow(H1b[s][kNext][w],2) + pow(H1c[s][kNext][w],2))\
                    + (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*
                    (pow(H2a[s][kNext][w],2) + pow(H2b[s][kNext][w],2) + pow(H2c[s][kNext][w],2))\
                    + (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*
                    (pow(H3a[s][kNext][w],2) + pow(H3b[s][kNext][w],2) + pow(H3c[s][kNext][w],2))\
                    + (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*
                    (pow(H4a[s][kNext][w],2) + pow(H4b[s][kNext][w],2) + pow(H4c[s][kNext][w],2))))/3.;


    DtH = Keff*(2*(2*(H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] +
                      H1c[s][k][w]*H2c[s][k][w])*
                   (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] +
                    H1c[s][k][wNext]*H2c[s][k][wNext]) +
                   2*(H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] +
                      H1c[s][k][w]*H3c[s][k][w])*
                   (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] +
                    H1c[s][k][wNext]*H3c[s][k][wNext]) +
                   2*(H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] +
                      H2c[s][k][w]*H3c[s][k][w])*
                   (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] +
                    H2c[s][k][wNext]*H3c[s][k][wNext]) +
                   2*(H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] +
                      H1c[s][k][w]*H4c[s][k][w])*
                   (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] +
                    H1c[s][k][wNext]*H4c[s][k][wNext]) +
                   2*(H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] +
                      H2c[s][k][w]*H4c[s][k][w])*
                   (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] +
                    H2c[s][k][wNext]*H4c[s][k][wNext]) +
                   2*(H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] +
                      H3c[s][k][w]*H4c[s][k][w])*
                   (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] +
                    H3c[s][k][wNext]*H4c[s][k][wNext]) +
                   (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*
                   (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2))\
                   + (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*
                   (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2))\
                   + (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*
                   (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2))\
                   + (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*
                   (pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2))))/3.;
    /*
    
    //---------
    DijHH = Geff*(2*(H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w]) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w]) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w]) +
                  2*(H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w]) +
                  2*(H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w]) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  2*(H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  2*(H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  2*(H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w]) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2)) +
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2)) +
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2))*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2))*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2))*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2))*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2))*(pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2))*(pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2)) +
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*(pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2)) +
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2))*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2))*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2))*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2))*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2))*(pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2))*(pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*(pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*(pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2))*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2))*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2))*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2))*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2))*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2))*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2))*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2))*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H1a[sNext][kPrev][w]*H2a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H2b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H2c[sNext][kPrev][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H1a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H2a[sNext][kPrev][w]*H3a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H3b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H3c[sNext][kPrev][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H1a[sNext][kPrev][w],2) + pow(H1b[sNext][kPrev][w],2) + pow(H1c[sNext][kPrev][w],2))*(pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H2a[sNext][kPrev][w],2) + pow(H2b[sNext][kPrev][w],2) + pow(H2c[sNext][kPrev][w],2))*(pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H3a[sNext][kPrev][w],2) + pow(H3b[sNext][kPrev][w],2) + pow(H3c[sNext][kPrev][w],2))*(pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H2a[sNext][k][w] + H1b[sNext][k][w]*H2b[sNext][k][w] + H1c[sNext][k][w]*H2c[sNext][k][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H3a[sNext][k][w] + H1b[sNext][k][w]*H3b[sNext][k][w] + H1c[sNext][k][w]*H3c[sNext][k][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H2a[sNext][k][w]*H3a[sNext][k][w] + H2b[sNext][k][w]*H3b[sNext][k][w] + H2c[sNext][k][w]*H3c[sNext][k][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[sNext][k][w],2) + pow(H1b[sNext][k][w],2) + pow(H1c[sNext][k][w],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[sNext][k][w],2) + pow(H2b[sNext][k][w],2) + pow(H2c[sNext][k][w],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[sNext][k][w],2) + pow(H3b[sNext][k][w],2) + pow(H3c[sNext][k][w],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) + 
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) + 
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) + 
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) + 
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*(pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) + 
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*(pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) + 
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2))*(pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) + 
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H1b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H1c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) + 
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H2b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H2c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) + 
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][kPrev][w]*H4a[sNext][kPrev][w] + H3b[sNext][kPrev][w]*H4b[sNext][kPrev][w] + H3c[sNext][kPrev][w]*H4c[sNext][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2)) + 
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[sNext][k][w]*H4a[sNext][k][w] + H1b[sNext][k][w]*H4b[sNext][k][w] + H1c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2))*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[sNext][k][w]*H4a[sNext][k][w] + H2b[sNext][k][w]*H4b[sNext][k][w] + H2c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2))*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[sNext][k][w]*H4a[sNext][k][w] + H3b[sNext][k][w]*H4b[sNext][k][w] + H3c[sNext][k][w]*H4c[sNext][k][w])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2))*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2))*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2))*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2))*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)) + 
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2))*
                  (pow(H4a[sNext][k][w],2) + pow(H4b[sNext][k][w],2) + pow(H4c[sNext][k][w],2))*
                  (pow(H4a[sNext][kPrev][w],2) + pow(H4b[sNext][kPrev][w],2) + pow(H4c[sNext][kPrev][w],2)));
    
    
    
    
    
    DtjHH = Geff*(2*(H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext]) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w]) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext]) +
                  2*(H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext]) +
                  2*(H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext]) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  2*(H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  2*(H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  2*(H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext]) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2)) +
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2)) +
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2))*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2))*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2))*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2))*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2))*(pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2)) +
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*(pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2)) +
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2))*(pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2))*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2))*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2))*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2))*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2))*(pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2))*(pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*(pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*(pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2))*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2))*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2))*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2))*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2))*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2))*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2))*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2))*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H1a[s][kPrev][wNext]*H2a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H2b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H2c[s][kPrev][wNext])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H1a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H2a[s][kPrev][wNext]*H3a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H3b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H3c[s][kPrev][wNext])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[s][kPrev][wNext],2) + pow(H1b[s][kPrev][wNext],2) + pow(H1c[s][kPrev][wNext],2))*(pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[s][kPrev][wNext],2) + pow(H2b[s][kPrev][wNext],2) + pow(H2c[s][kPrev][wNext],2))*(pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[s][kPrev][wNext],2) + pow(H3b[s][kPrev][wNext],2) + pow(H3c[s][kPrev][wNext],2))*(pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H1a[s][kPrev][w]*H2a[s][kPrev][w] + H1b[s][kPrev][w]*H2b[s][kPrev][w] + H1c[s][kPrev][w]*H2c[s][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H1a[s][kPrev][w]*H3a[s][kPrev][w] + H1b[s][kPrev][w]*H3b[s][kPrev][w] + H1c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H2a[s][kPrev][w]*H3a[s][kPrev][w] + H2b[s][kPrev][w]*H3b[s][kPrev][w] + H2c[s][kPrev][w]*H3c[s][kPrev][w])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][kPrev][w],2) + pow(H1b[s][kPrev][w],2) + pow(H1c[s][kPrev][w],2))*(pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][kPrev][w],2) + pow(H2b[s][kPrev][w],2) + pow(H2c[s][kPrev][w],2))*(pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][kPrev][w],2) + pow(H3b[s][kPrev][w],2) + pow(H3c[s][kPrev][w],2))*(pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) +
                  (H1a[s][k][wNext]*H2a[s][k][wNext] + H1b[s][k][wNext]*H2b[s][k][wNext] + H1c[s][k][wNext]*H2c[s][k][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) + 
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) + 
                  (H1a[s][k][wNext]*H3a[s][k][wNext] + H1b[s][k][wNext]*H3b[s][k][wNext] + H1c[s][k][wNext]*H3c[s][k][wNext])*
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) + 
                  (H2a[s][k][wNext]*H3a[s][k][wNext] + H2b[s][k][wNext]*H3b[s][k][wNext] + H2c[s][k][wNext]*H3c[s][k][wNext])*
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) + 
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H1a[s][k][wNext],2) + pow(H1b[s][k][wNext],2) + pow(H1c[s][k][wNext],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) + 
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H2a[s][k][wNext],2) + pow(H2b[s][k][wNext],2) + pow(H2c[s][k][wNext],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) + 
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H3a[s][k][wNext],2) + pow(H3b[s][k][wNext],2) + pow(H3c[s][k][wNext],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) + 
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H1b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H1c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) + 
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H2b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H2c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) + 
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][wNext]*H4a[s][kPrev][wNext] + H3b[s][kPrev][wNext]*H4b[s][kPrev][wNext] + H3c[s][kPrev][wNext]*H4c[s][kPrev][wNext])*
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2)) + 
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H1a[s][k][w]*H2a[s][k][w] + H1b[s][k][w]*H2b[s][k][w] + H1c[s][k][w]*H2c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H1a[s][k][w]*H3a[s][k][w] + H1b[s][k][w]*H3b[s][k][w] + H1c[s][k][w]*H3c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H2a[s][k][w]*H3a[s][k][w] + H2b[s][k][w]*H3b[s][k][w] + H2c[s][k][w]*H3c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H1a[s][k][w],2) + pow(H1b[s][k][w],2) + pow(H1c[s][k][w],2))*(pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H2a[s][k][w],2) + pow(H2b[s][k][w],2) + pow(H2c[s][k][w],2))*(pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H3a[s][k][w],2) + pow(H3b[s][k][w],2) + pow(H3c[s][k][w],2))*(pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][kPrev][w]*H4a[s][kPrev][w] + H1b[s][kPrev][w]*H4b[s][kPrev][w] + H1c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2))*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][kPrev][w]*H4a[s][kPrev][w] + H2b[s][kPrev][w]*H4b[s][kPrev][w] + H2c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2))*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][kPrev][w]*H4a[s][kPrev][w] + H3b[s][kPrev][w]*H4b[s][kPrev][w] + H3c[s][kPrev][w]*H4c[s][kPrev][w])*
                  (pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2))*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H1a[s][k][w]*H4a[s][k][w] + H1b[s][k][w]*H4b[s][k][w] + H1c[s][k][w]*H4c[s][k][w])*
                  (H1a[s][k][wNext]*H4a[s][k][wNext] + H1b[s][k][wNext]*H4b[s][k][wNext] + H1c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2))*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H2a[s][k][w]*H4a[s][k][w] + H2b[s][k][w]*H4b[s][k][w] + H2c[s][k][w]*H4c[s][k][w])*
                  (H2a[s][k][wNext]*H4a[s][k][wNext] + H2b[s][k][wNext]*H4b[s][k][wNext] + H2c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2))*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (H3a[s][k][w]*H4a[s][k][w] + H3b[s][k][w]*H4b[s][k][w] + H3c[s][k][w]*H4c[s][k][w])*
                  (H3a[s][k][wNext]*H4a[s][k][wNext] + H3b[s][k][wNext]*H4b[s][k][wNext] + H3c[s][k][wNext]*H4c[s][k][wNext])*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2))*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)) + 
                  (pow(H4a[s][k][w],2) + pow(H4b[s][k][w],2) + pow(H4c[s][k][w],2))*(pow(H4a[s][k][wNext],2) + pow(H4b[s][k][wNext],2) + pow(H4c[s][k][wNext],2))*
                  (pow(H4a[s][kPrev][w],2) + pow(H4b[s][kPrev][w],2) + pow(H4c[s][kPrev][w],2))*
                  (pow(H4a[s][kPrev][wNext],2) + pow(H4b[s][kPrev][wNext],2) + pow(H4c[s][kPrev][wNext],2)));

    
*/
    
    eConfig = VH + DiH + DjH  + DtH; // + DijHH +DtjHH;
    
    return eConfig;
}



std::complex<double> E0(int i, int j, int t){
    std::complex<double> e0=0. + 0.*I1;
    

    H1a0[i][j][t]=H1a[i][j][t];
    H1b0[i][j][t]=H1b[i][j][t];
    H1c0[i][j][t]=H1c[i][j][t];
    
    H2a0[i][j][t]=H2a[i][j][t];
    H2b0[i][j][t]=H2b[i][j][t];
    H2c0[i][j][t]=H2c[i][j][t];
    
    H3a0[i][j][t]=H3a[i][j][t];
    H3b0[i][j][t]=H3b[i][j][t];
    H3c0[i][j][t]=H3c[i][j][t];
    
    H4a0[i][j][t]=H4a[i][j][t];
    H4b0[i][j][t]=H4b[i][j][t];
    H4c0[i][j][t]=H4c[i][j][t];
    
    
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
            
            int s = IarrayNum[sindex];
            int k = JarrayNum[kindex];
            int w = TarrayNum[windex];
            
            int sPrev = s == 0 ? Lx-1 : s-1;
            int sNext = s == Lx-1 ? 0 : s+1;
            int kPrev = k == 0 ? Ly-1 : k-1;
            int kNext = k == Ly-1 ? 0 : k+1;
            int wPrev = w == 0 ? Lt-1 : w-1;
            int wNext = w == Lt-1 ? 0 : w+1;
    

            

 
            e0 += EConfiguration(s,k,w);
            

             
            }}}

    
    return e0;}






std::complex<double> EfH1(int i, int j, int t){
    std::complex<double> eflip1=0. + 0.*I1;

    
    a=rand()%2;
    
    // Update Higgs Field Amplitude
    
    double normH1;
    
    normH1= sqrt(H1a[i][j][t]*H1a[i][j][t] + H1b[i][j][t]*H1b[i][j][t] + H1c[i][j][t]*H1c[i][j][t]);
    AmpH1 = (normH1 + d_Amp*(1.-2.*std_rand()))/normH1;     // generate random amplitudes for Higgs H1
    
    
    
    H1aflip[i][j][t] =AmpH1*H1a[i][j][t];
    H1bflip[i][j][t] =AmpH1*H1b[i][j][t];
    H1cflip[i][j][t] =AmpH1*H1c[i][j][t];
    
    
    
    
    
    
    H1a[i][j][t]=H1aflip[i][j][t];
    H1b[i][j][t]=H1bflip[i][j][t];
    H1c[i][j][t]=H1cflip[i][j][t];
    
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H1a[s][k][w]=H1aflip[s][k][w];
                    H1b[s][k][w]=H1bflip[s][k][w];
                    H1c[s][k][w]=H1cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H1a[sPrev][k][w]=H1aflip[sPrev][k][w];
                    H1b[sPrev][k][w]=H1bflip[sPrev][k][w];
                    H1c[sPrev][k][w]=H1cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H1a[s][kPrev][w]=H1aflip[s][kPrev][w];
                    H1b[s][kPrev][w]=H1bflip[s][kPrev][w];
                    H1c[s][kPrev][w]=H1cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H1a[sNext][k][w]=H1aflip[sNext][k][w];
                    H1b[sNext][k][w]=H1bflip[sNext][k][w];
                    H1c[sNext][k][w]=H1cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H1a[s][kNext][w]=H1aflip[s][kNext][w];
                    H1b[s][kNext][w]=H1bflip[s][kNext][w];
                    H1c[s][kNext][w]=H1cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H1a[s][k][wPrev]=H1aflip[s][k][wPrev];
                    H1b[s][k][wPrev]=H1bflip[s][k][wPrev];
                    H1c[s][k][wPrev]=H1cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H1a[s][k][wNext]=H1aflip[s][k][wNext];
                    H1b[s][k][wNext]=H1bflip[s][k][wNext];
                    H1c[s][k][wNext]=H1cflip[s][k][wNext];}
                
                
                
            
            
            eflip1 += EConfiguration(s,k,w);
    
            }}}
    

    
    return eflip1;}


std::complex<double> EfH2(int i, int j, int t){
    std::complex<double> eflip2=0. + 0.*I1;
    
    
    a=rand()%2;
    
    // Update Higgs Field Amplitude
    
    double normH2;
    
    normH2= sqrt(H2a[i][j][t]*H2a[i][j][t] + H2b[i][j][t]*H2b[i][j][t] + H2c[i][j][t]*H2c[i][j][t]);
    AmpH2 = (normH2 + d_Amp*(1.-2.*std_rand()))/normH2;     // generate random amplitudes for Higgs H2
    
    
    
    H2aflip[i][j][t] =AmpH2*H2a[i][j][t];
    H2bflip[i][j][t] =AmpH2*H2b[i][j][t];
    H2cflip[i][j][t] =AmpH2*H2c[i][j][t];
    
    
    
    
    
    
    H2a[i][j][t]=H2aflip[i][j][t];
    H2b[i][j][t]=H2bflip[i][j][t];
    H2c[i][j][t]=H2cflip[i][j][t];
    
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H2a[s][k][w]=H2aflip[s][k][w];
                    H2b[s][k][w]=H2bflip[s][k][w];
                    H2c[s][k][w]=H2cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H2a[sPrev][k][w]=H2aflip[sPrev][k][w];
                    H2b[sPrev][k][w]=H2bflip[sPrev][k][w];
                    H2c[sPrev][k][w]=H2cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H2a[s][kPrev][w]=H2aflip[s][kPrev][w];
                    H2b[s][kPrev][w]=H2bflip[s][kPrev][w];
                    H2c[s][kPrev][w]=H2cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H2a[sNext][k][w]=H2aflip[sNext][k][w];
                    H2b[sNext][k][w]=H2bflip[sNext][k][w];
                    H2c[sNext][k][w]=H2cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H2a[s][kNext][w]=H2aflip[s][kNext][w];
                    H2b[s][kNext][w]=H2bflip[s][kNext][w];
                    H2c[s][kNext][w]=H2cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H2a[s][k][wPrev]=H2aflip[s][k][wPrev];
                    H2b[s][k][wPrev]=H2bflip[s][k][wPrev];
                    H2c[s][k][wPrev]=H2cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H2a[s][k][wNext]=H2aflip[s][k][wNext];
                    H2b[s][k][wNext]=H2bflip[s][k][wNext];
                    H2c[s][k][wNext]=H2cflip[s][k][wNext];}
                
                
                
                
                //Higher order derivatives
                
//                else if(sPrev==i && kPrev==j && w==t){
//                    H2a[sPrev][kPrev][w]=H2aflip[sPrev][kPrev][w];
//                    H2b[sPrev][kPrev][w]=H2bflip[sPrev][kPrev][w];
//                    H2c[sPrev][kPrev][w]=H2cflip[sPrev][kPrev][w];}
//                else if(sNext==i && kNext==j && w==t){
//                    H2a[sNext][kNext][w]=H2aflip[sNext][kNext][w];
//                    H2b[sNext][kNext][w]=H2bflip[sNext][kNext][w];
//                    H2c[sNext][kNext][w]=H2cflip[sNext][kNext][w];}
//                
//                else if(sNext==i && kPrev==j && w==t){
//                    H2a[sNext][kPrev][w]=H2aflip[sNext][kPrev][w];
//                    H2b[sNext][kPrev][w]=H2bflip[sNext][kPrev][w];
//                    H2c[sNext][kPrev][w]=H2cflip[sNext][kPrev][w];}
//                
//                else if(sPrev==i && kNext==j && w==t){
//                    H2a[sPrev][kNext][w]=H2aflip[sPrev][kNext][w];
//                    H2b[sPrev][kNext][w]=H2bflip[sPrev][kNext][w];
//                    H2c[sPrev][kNext][w]=H2cflip[sPrev][kNext][w];}
//                else if(sPrev==i && k==j && wPrev==t){
//                    H2a[sPrev][k][wPrev]=H2aflip[sPrev][k][wPrev];
//                    H2b[sPrev][k][wPrev]=H2bflip[sPrev][k][wPrev];
//                    H2c[sPrev][k][wPrev]=H2cflip[sPrev][k][wPrev];}
//                
//                else if(sPrev==i && k==j && wNext==t){
//                    H2a[sPrev][k][wNext]=H2aflip[sPrev][k][wNext];
//                    H2b[sPrev][k][wNext]=H2bflip[sPrev][k][wNext];
//                    H2c[sPrev][k][wNext]=H2cflip[sPrev][k][wNext];}
//                
//                else if(sNext==i && k==j && wPrev==t){
//                    H2a[sNext][k][wPrev]=H2aflip[sNext][k][wPrev];
//                    H2b[sNext][k][wPrev]=H2bflip[sNext][k][wPrev];
//                    H2c[sNext][k][wPrev]=H2cflip[sNext][k][wPrev];}
//                
//                else if(sNext==i && k==j && wNext==t){
//                    H2a[sNext][k][wNext]=H2aflip[sNext][k][wNext];
//                    H2b[sNext][k][wNext]=H2bflip[sNext][k][wNext];
//                    H2c[sNext][k][wNext]=H2cflip[sNext][k][wNext];}
//                
//                
//                
//                else if(s==i && kPrev==j && wPrev==t){
//                    H2a[s][kPrev][wPrev]=H2aflip[s][kPrev][wPrev];
//                    H2b[s][kPrev][wPrev]=H2bflip[s][kPrev][wPrev];
//                    H2c[s][kPrev][wPrev]=H2cflip[s][kPrev][wPrev];}
//                
//                else if(s==i && kPrev==j && wNext==t){
//                    H2a[s][kPrev][wNext]=H2aflip[s][kPrev][wNext];
//                    H2b[s][kPrev][wNext]=H2bflip[s][kPrev][wNext];
//                    H2c[s][kPrev][wNext]=H2cflip[s][kPrev][wNext];}
//                
//                
//                else if(s==i && kNext==j && wPrev==t){
//                    H2a[s][kNext][wPrev]=H2aflip[s][kNext][wPrev];
//                    H2b[s][kNext][wPrev]=H2bflip[s][kNext][wPrev];
//                    H2c[s][kNext][wPrev]=H2cflip[s][kNext][wPrev];}
//                
//                else if(s==i && kNext==j && wNext==t){
//                    H2a[s][kNext][wNext]=H2aflip[s][kNext][wNext];
//                    H2b[s][kNext][wNext]=H2bflip[s][kNext][wNext];
//                    H2c[s][kNext][wNext]=H2cflip[s][kNext][wNext];}
                

            
            
            
            eflip2+= EConfiguration(s,k,w);
            }}}

    
    return eflip2;}


std::complex<double> EfH3(int i, int j, int t){
    std::complex<double> eflip3=0. + 0.*I1;
    

    a=rand()%2;
    
    // Update Higgs Field Amplitude
    
    double normH3;
    
    normH3= sqrt(H3a[i][j][t]*H3a[i][j][t] + H3b[i][j][t]*H3b[i][j][t] + H3c[i][j][t]*H3c[i][j][t]);
    AmpH3 = (normH3 + d_Amp*(1.-2.*std_rand()))/normH3;     // generate random amplitudes for Higgs H3
    
    
    
    H3aflip[i][j][t] =AmpH3*H3a[i][j][t];
    H3bflip[i][j][t] =AmpH3*H3b[i][j][t];
    H3cflip[i][j][t] =AmpH3*H3c[i][j][t];
    
    
    
    
    
    
    H3a[i][j][t]=H3aflip[i][j][t];
    H3b[i][j][t]=H3bflip[i][j][t];
    H3c[i][j][t]=H3cflip[i][j][t];
    
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H3a[s][k][w]=H3aflip[s][k][w];
                    H3b[s][k][w]=H3bflip[s][k][w];
                    H3c[s][k][w]=H3cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H3a[sPrev][k][w]=H3aflip[sPrev][k][w];
                    H3b[sPrev][k][w]=H3bflip[sPrev][k][w];
                    H3c[sPrev][k][w]=H3cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H3a[s][kPrev][w]=H3aflip[s][kPrev][w];
                    H3b[s][kPrev][w]=H3bflip[s][kPrev][w];
                    H3c[s][kPrev][w]=H3cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H3a[sNext][k][w]=H3aflip[sNext][k][w];
                    H3b[sNext][k][w]=H3bflip[sNext][k][w];
                    H3c[sNext][k][w]=H3cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H3a[s][kNext][w]=H3aflip[s][kNext][w];
                    H3b[s][kNext][w]=H3bflip[s][kNext][w];
                    H3c[s][kNext][w]=H3cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H3a[s][k][wPrev]=H3aflip[s][k][wPrev];
                    H3b[s][k][wPrev]=H3bflip[s][k][wPrev];
                    H3c[s][k][wPrev]=H3cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H3a[s][k][wNext]=H3aflip[s][k][wNext];
                    H3b[s][k][wNext]=H3bflip[s][k][wNext];
                    H3c[s][k][wNext]=H3cflip[s][k][wNext];}
                
                

    
    
    
            
            
            eflip3+=EConfiguration(s,k,w);
            }}}

    
    return eflip3;}


std::complex<double> EfH4(int i, int j, int t){
    std::complex<double> eflip4=0. + 0.*I1;
    
    
    a=rand()%2;
    
    // Update Higgs Field Amplitude
    
    double normH4;
    
    normH4= sqrt(H4a[i][j][t]*H4a[i][j][t] + H4b[i][j][t]*H4b[i][j][t] + H4c[i][j][t]*H4c[i][j][t]);
    AmpH4 = (normH4 + d_Amp*(1.-2.*std_rand()))/normH4;     // generate random amplitudes for Higgs H4
    
    
    
    H4aflip[i][j][t] =AmpH4*H4a[i][j][t];
    H4bflip[i][j][t] =AmpH4*H4b[i][j][t];
    H4cflip[i][j][t] =AmpH4*H4c[i][j][t];
    
    
    
    
    
    
    H4a[i][j][t]=H4aflip[i][j][t];
    H4b[i][j][t]=H4bflip[i][j][t];
    H4c[i][j][t]=H4cflip[i][j][t];
    
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H4a[s][k][w]=H4aflip[s][k][w];
                    H4b[s][k][w]=H4bflip[s][k][w];
                    H4c[s][k][w]=H4cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H4a[sPrev][k][w]=H4aflip[sPrev][k][w];
                    H4b[sPrev][k][w]=H4bflip[sPrev][k][w];
                    H4c[sPrev][k][w]=H4cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H4a[s][kPrev][w]=H4aflip[s][kPrev][w];
                    H4b[s][kPrev][w]=H4bflip[s][kPrev][w];
                    H4c[s][kPrev][w]=H4cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H4a[sNext][k][w]=H4aflip[sNext][k][w];
                    H4b[sNext][k][w]=H4bflip[sNext][k][w];
                    H4c[sNext][k][w]=H4cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H4a[s][kNext][w]=H4aflip[s][kNext][w];
                    H4b[s][kNext][w]=H4bflip[s][kNext][w];
                    H4c[s][kNext][w]=H4cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H4a[s][k][wPrev]=H4aflip[s][k][wPrev];
                    H4b[s][k][wPrev]=H4bflip[s][k][wPrev];
                    H4c[s][k][wPrev]=H4cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H4a[s][k][wNext]=H4aflip[s][k][wNext];
                    H4b[s][k][wNext]=H4bflip[s][k][wNext];
                    H4c[s][k][wNext]=H4cflip[s][k][wNext];}
                
                
                
            
            
            eflip4+=EConfiguration(s,k,w);
            }}}
    

    
    return eflip4;}






std::complex<double> EfH1Rab(int i, int j, int t){
    std::complex<double> eflip1=0. + 0.*I1;
    
    
    // Update Higgs Fields via small angle rotations
    
    
    alpha = (1.-2.*std_rand())*d_Phi;
    C_alpha = cos(alpha);            // here alpha_max = PI
    S_alpha = sin(alpha);
    
    
    
    H1aflip[i][j][t] = C_alpha*H1a[i][j][t] + S_alpha*H1b[i][j][t];
    H1bflip[i][j][t] = -S_alpha*H1a[i][j][t] + C_alpha*H1b[i][j][t];
    H1cflip[i][j][t] =H1c[i][j][t];
    
    
    
    H1a[i][j][t]=H1aflip[i][j][t];
    H1b[i][j][t]=H1bflip[i][j][t];
    H1c[i][j][t]=H1cflip[i][j][t];
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H1a[s][k][w]=H1aflip[s][k][w];
                    H1b[s][k][w]=H1bflip[s][k][w];
                    H1c[s][k][w]=H1cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H1a[sPrev][k][w]=H1aflip[sPrev][k][w];
                    H1b[sPrev][k][w]=H1bflip[sPrev][k][w];
                    H1c[sPrev][k][w]=H1cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H1a[s][kPrev][w]=H1aflip[s][kPrev][w];
                    H1b[s][kPrev][w]=H1bflip[s][kPrev][w];
                    H1c[s][kPrev][w]=H1cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H1a[sNext][k][w]=H1aflip[sNext][k][w];
                    H1b[sNext][k][w]=H1bflip[sNext][k][w];
                    H1c[sNext][k][w]=H1cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H1a[s][kNext][w]=H1aflip[s][kNext][w];
                    H1b[s][kNext][w]=H1bflip[s][kNext][w];
                    H1c[s][kNext][w]=H1cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H1a[s][k][wPrev]=H1aflip[s][k][wPrev];
                    H1b[s][k][wPrev]=H1bflip[s][k][wPrev];
                    H1c[s][k][wPrev]=H1cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H1a[s][k][wNext]=H1aflip[s][k][wNext];
                    H1b[s][k][wNext]=H1bflip[s][k][wNext];
                    H1c[s][k][wNext]=H1cflip[s][k][wNext];}
                
                
                
            
            
            eflip1+= EConfiguration(s,k,w);
            }}}
    
    return eflip1;}

std::complex<double> EfH1Rbc(int i, int j, int t){
    std::complex<double> eflip1=0. + 0.*I1;
    
    
    // Update Higgs Fields via small angle rotations
    
    
    alpha = (1.-2.*std_rand())*d_Phi;
    C_alpha = cos(alpha);            // here alpha_max = PI
    S_alpha = sin(alpha);
    
    
    
    H1bflip[i][j][t] = C_alpha*H1b[i][j][t] + S_alpha*H1c[i][j][t];
    H1cflip[i][j][t] = -S_alpha*H1b[i][j][t] + C_alpha*H1c[i][j][t];
    H1aflip[i][j][t] =H1a[i][j][t];
    
    
    H1a[i][j][t]=H1aflip[i][j][t];
    H1b[i][j][t]=H1bflip[i][j][t];
    H1c[i][j][t]=H1cflip[i][j][t];
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H1a[s][k][w]=H1aflip[s][k][w];
                    H1b[s][k][w]=H1bflip[s][k][w];
                    H1c[s][k][w]=H1cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H1a[sPrev][k][w]=H1aflip[sPrev][k][w];
                    H1b[sPrev][k][w]=H1bflip[sPrev][k][w];
                    H1c[sPrev][k][w]=H1cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H1a[s][kPrev][w]=H1aflip[s][kPrev][w];
                    H1b[s][kPrev][w]=H1bflip[s][kPrev][w];
                    H1c[s][kPrev][w]=H1cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H1a[sNext][k][w]=H1aflip[sNext][k][w];
                    H1b[sNext][k][w]=H1bflip[sNext][k][w];
                    H1c[sNext][k][w]=H1cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H1a[s][kNext][w]=H1aflip[s][kNext][w];
                    H1b[s][kNext][w]=H1bflip[s][kNext][w];
                    H1c[s][kNext][w]=H1cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H1a[s][k][wPrev]=H1aflip[s][k][wPrev];
                    H1b[s][k][wPrev]=H1bflip[s][k][wPrev];
                    H1c[s][k][wPrev]=H1cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H1a[s][k][wNext]=H1aflip[s][k][wNext];
                    H1b[s][k][wNext]=H1bflip[s][k][wNext];
                    H1c[s][k][wNext]=H1cflip[s][k][wNext];}
                
                
                
            
            
            eflip1+= EConfiguration(s,k,w);
            }}}
    
    return eflip1;}





std::complex<double> EfH2Rab(int i, int j, int t){
    std::complex<double> eflip1=0. + 0.*I1;
    
    
    // Update Higgs Fields via small angle rotations
    
    
    alpha = (1.-2.*std_rand())*d_Phi;
    C_alpha = cos(alpha);            // here alpha_max = PI
    S_alpha = sin(alpha);
    
    
    
    H2aflip[i][j][t] = C_alpha*H2a[i][j][t] + S_alpha*H2b[i][j][t];
    H2bflip[i][j][t] = -S_alpha*H2a[i][j][t] + C_alpha*H2b[i][j][t];
    H2cflip[i][j][t] =H2c[i][j][t];
    
    
    
    H2a[i][j][t]=H2aflip[i][j][t];
    H2b[i][j][t]=H2bflip[i][j][t];
    H2c[i][j][t]=H2cflip[i][j][t];
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H2a[s][k][w]=H2aflip[s][k][w];
                    H2b[s][k][w]=H2bflip[s][k][w];
                    H2c[s][k][w]=H2cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H2a[sPrev][k][w]=H2aflip[sPrev][k][w];
                    H2b[sPrev][k][w]=H2bflip[sPrev][k][w];
                    H2c[sPrev][k][w]=H2cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H2a[s][kPrev][w]=H2aflip[s][kPrev][w];
                    H2b[s][kPrev][w]=H2bflip[s][kPrev][w];
                    H2c[s][kPrev][w]=H2cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H2a[sNext][k][w]=H2aflip[sNext][k][w];
                    H2b[sNext][k][w]=H2bflip[sNext][k][w];
                    H2c[sNext][k][w]=H2cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H2a[s][kNext][w]=H2aflip[s][kNext][w];
                    H2b[s][kNext][w]=H2bflip[s][kNext][w];
                    H2c[s][kNext][w]=H2cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H2a[s][k][wPrev]=H2aflip[s][k][wPrev];
                    H2b[s][k][wPrev]=H2bflip[s][k][wPrev];
                    H2c[s][k][wPrev]=H2cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H2a[s][k][wNext]=H2aflip[s][k][wNext];
                    H2b[s][k][wNext]=H2bflip[s][k][wNext];
                    H2c[s][k][wNext]=H2cflip[s][k][wNext];}
                
                
                
            
            
            eflip1+= EConfiguration(s,k,w);
            }}}
    
    
    
    return eflip1;}

std::complex<double> EfH2Rbc(int i, int j, int t){
    std::complex<double> eflip1=0. + 0.*I1;
    
    
    // Update Higgs Fields via small angle rotations
    
    
    alpha = (1.-2.*std_rand())*d_Phi;
    C_alpha = cos(alpha);            // here alpha_max = PI
    S_alpha = sin(alpha);
    
    
    
    H2bflip[i][j][t] = C_alpha*H2b[i][j][t] + S_alpha*H2c[i][j][t];
    H2cflip[i][j][t] = -S_alpha*H2b[i][j][t] + C_alpha*H2c[i][j][t];
    H2aflip[i][j][t] =H2a[i][j][t];
    
    
    H2a[i][j][t]=H2aflip[i][j][t];
    H2b[i][j][t]=H2bflip[i][j][t];
    H2c[i][j][t]=H2cflip[i][j][t];
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H2a[s][k][w]=H2aflip[s][k][w];
                    H2b[s][k][w]=H2bflip[s][k][w];
                    H2c[s][k][w]=H2cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H2a[sPrev][k][w]=H2aflip[sPrev][k][w];
                    H2b[sPrev][k][w]=H2bflip[sPrev][k][w];
                    H2c[sPrev][k][w]=H2cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H2a[s][kPrev][w]=H2aflip[s][kPrev][w];
                    H2b[s][kPrev][w]=H2bflip[s][kPrev][w];
                    H2c[s][kPrev][w]=H2cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H2a[sNext][k][w]=H2aflip[sNext][k][w];
                    H2b[sNext][k][w]=H2bflip[sNext][k][w];
                    H2c[sNext][k][w]=H2cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H2a[s][kNext][w]=H2aflip[s][kNext][w];
                    H2b[s][kNext][w]=H2bflip[s][kNext][w];
                    H2c[s][kNext][w]=H2cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H2a[s][k][wPrev]=H2aflip[s][k][wPrev];
                    H2b[s][k][wPrev]=H2bflip[s][k][wPrev];
                    H2c[s][k][wPrev]=H2cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H2a[s][k][wNext]=H2aflip[s][k][wNext];
                    H2b[s][k][wNext]=H2bflip[s][k][wNext];
                    H2c[s][k][wNext]=H2cflip[s][k][wNext];}
                
                
                
            
            
            
            
            eflip1+= EConfiguration(s,k,w);
            }}}
    
    
    
    return eflip1;}





std::complex<double> EfH3Rab(int i, int j, int t){
    std::complex<double> eflip1=0. + 0.*I1;
    
    
    // Update Higgs Fields via small angle rotations
    
    
    alpha = (1.-2.*std_rand())*d_Phi;
    C_alpha = cos(alpha);            // here alpha_max = PI
    S_alpha = sin(alpha);
    
    
    
    H3aflip[i][j][t] = C_alpha*H3a[i][j][t] + S_alpha*H3b[i][j][t];
    H3bflip[i][j][t] = -S_alpha*H3a[i][j][t] + C_alpha*H3b[i][j][t];
    H3cflip[i][j][t] =H3c[i][j][t];
    
    
    
    H3a[i][j][t]=H3aflip[i][j][t];
    H3b[i][j][t]=H3bflip[i][j][t];
    H3c[i][j][t]=H3cflip[i][j][t];
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H3a[s][k][w]=H3aflip[s][k][w];
                    H3b[s][k][w]=H3bflip[s][k][w];
                    H3c[s][k][w]=H3cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H3a[sPrev][k][w]=H3aflip[sPrev][k][w];
                    H3b[sPrev][k][w]=H3bflip[sPrev][k][w];
                    H3c[sPrev][k][w]=H3cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H3a[s][kPrev][w]=H3aflip[s][kPrev][w];
                    H3b[s][kPrev][w]=H3bflip[s][kPrev][w];
                    H3c[s][kPrev][w]=H3cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H3a[sNext][k][w]=H3aflip[sNext][k][w];
                    H3b[sNext][k][w]=H3bflip[sNext][k][w];
                    H3c[sNext][k][w]=H3cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H3a[s][kNext][w]=H3aflip[s][kNext][w];
                    H3b[s][kNext][w]=H3bflip[s][kNext][w];
                    H3c[s][kNext][w]=H3cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H3a[s][k][wPrev]=H3aflip[s][k][wPrev];
                    H3b[s][k][wPrev]=H3bflip[s][k][wPrev];
                    H3c[s][k][wPrev]=H3cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H3a[s][k][wNext]=H3aflip[s][k][wNext];
                    H3b[s][k][wNext]=H3bflip[s][k][wNext];
                    H3c[s][k][wNext]=H3cflip[s][k][wNext];}
                
                
                
            
            
            eflip1+= EConfiguration(s,k,w);
            }}}
    
    
    
    return eflip1;}

std::complex<double> EfH3Rbc(int i, int j, int t){
    std::complex<double> eflip1=0. + 0.*I1;
    
    
    // Update Higgs Fields via small angle rotations
    
    
    alpha = (1.-2.*std_rand())*d_Phi;
    C_alpha = cos(alpha);            // here alpha_max = PI
    S_alpha = sin(alpha);
    
    
    
    H3bflip[i][j][t] = C_alpha*H3b[i][j][t] + S_alpha*H3c[i][j][t];
    H3cflip[i][j][t] = -S_alpha*H3b[i][j][t] + C_alpha*H3c[i][j][t];
    H3aflip[i][j][t] =H3a[i][j][t];
    
    
    H3a[i][j][t]=H3aflip[i][j][t];
    H3b[i][j][t]=H3bflip[i][j][t];
    H3c[i][j][t]=H3cflip[i][j][t];
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H3a[s][k][w]=H3aflip[s][k][w];
                    H3b[s][k][w]=H3bflip[s][k][w];
                    H3c[s][k][w]=H3cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H3a[sPrev][k][w]=H3aflip[sPrev][k][w];
                    H3b[sPrev][k][w]=H3bflip[sPrev][k][w];
                    H3c[sPrev][k][w]=H3cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H3a[s][kPrev][w]=H3aflip[s][kPrev][w];
                    H3b[s][kPrev][w]=H3bflip[s][kPrev][w];
                    H3c[s][kPrev][w]=H3cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H3a[sNext][k][w]=H3aflip[sNext][k][w];
                    H3b[sNext][k][w]=H3bflip[sNext][k][w];
                    H3c[sNext][k][w]=H3cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H3a[s][kNext][w]=H3aflip[s][kNext][w];
                    H3b[s][kNext][w]=H3bflip[s][kNext][w];
                    H3c[s][kNext][w]=H3cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H3a[s][k][wPrev]=H3aflip[s][k][wPrev];
                    H3b[s][k][wPrev]=H3bflip[s][k][wPrev];
                    H3c[s][k][wPrev]=H3cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H3a[s][k][wNext]=H3aflip[s][k][wNext];
                    H3b[s][k][wNext]=H3bflip[s][k][wNext];
                    H3c[s][k][wNext]=H3cflip[s][k][wNext];}
                
                
                

            
            
            
            
            
            eflip1+= EConfiguration(s,k,w);
            }}}
    
    
    
    return eflip1;}





std::complex<double> EfH4Rab(int i, int j, int t){
    std::complex<double> eflip1=0. + 0.*I1;
    
    
    // Update Higgs Fields via small angle rotations
    
    
    alpha = (1.-2.*std_rand())*d_Phi;
    C_alpha = cos(alpha);            // here alpha_max = PI
    S_alpha = sin(alpha);
    
    
    
    H4aflip[i][j][t] = C_alpha*H4a[i][j][t] + S_alpha*H4b[i][j][t];
    H4bflip[i][j][t] = -S_alpha*H4a[i][j][t] + C_alpha*H4b[i][j][t];
    H4cflip[i][j][t] =H4c[i][j][t];
    
    
    
    H4a[i][j][t]=H4aflip[i][j][t];
    H4b[i][j][t]=H4bflip[i][j][t];
    H4c[i][j][t]=H4cflip[i][j][t];
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H4a[s][k][w]=H4aflip[s][k][w];
                    H4b[s][k][w]=H4bflip[s][k][w];
                    H4c[s][k][w]=H4cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H4a[sPrev][k][w]=H4aflip[sPrev][k][w];
                    H4b[sPrev][k][w]=H4bflip[sPrev][k][w];
                    H4c[sPrev][k][w]=H4cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H4a[s][kPrev][w]=H4aflip[s][kPrev][w];
                    H4b[s][kPrev][w]=H4bflip[s][kPrev][w];
                    H4c[s][kPrev][w]=H4cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H4a[sNext][k][w]=H4aflip[sNext][k][w];
                    H4b[sNext][k][w]=H4bflip[sNext][k][w];
                    H4c[sNext][k][w]=H4cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H4a[s][kNext][w]=H4aflip[s][kNext][w];
                    H4b[s][kNext][w]=H4bflip[s][kNext][w];
                    H4c[s][kNext][w]=H4cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H4a[s][k][wPrev]=H4aflip[s][k][wPrev];
                    H4b[s][k][wPrev]=H4bflip[s][k][wPrev];
                    H4c[s][k][wPrev]=H4cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H4a[s][k][wNext]=H4aflip[s][k][wNext];
                    H4b[s][k][wNext]=H4bflip[s][k][wNext];
                    H4c[s][k][wNext]=H4cflip[s][k][wNext];}
                
                
                
            
            
            eflip1+= EConfiguration(s,k,w);
            }}}
    
    
    
    return eflip1;}

std::complex<double> EfH4Rbc(int i, int j, int t){
    std::complex<double> eflip1=0. + 0.*I1;
    
    
    // Update Higgs Fields via small angle rotations
    
    
    alpha = (1.-2.*std_rand())*d_Phi;
    C_alpha = cos(alpha);            // here alpha_max = PI
    S_alpha = sin(alpha);
    
    
    
    H4bflip[i][j][t] = C_alpha*H4b[i][j][t] + S_alpha*H4c[i][j][t];
    H4cflip[i][j][t] = -S_alpha*H4b[i][j][t] + C_alpha*H4c[i][j][t];
    H4aflip[i][j][t] =H4a[i][j][t];
    
    
    H4a[i][j][t]=H4aflip[i][j][t];
    H4b[i][j][t]=H4bflip[i][j][t];
    H4c[i][j][t]=H4cflip[i][j][t];
    
    int iPrev = i == 0 ? Lx-1 : i-1;
    int iNext = i == Lx-1 ? 0 : i+1;
    int jPrev = j == 0 ? Ly-1 : j-1;
    int jNext = j == Ly-1 ? 0 : j+1;
    int tPrev = t == 0 ? Lt-1 : t-1;
    int tNext = t == Lt-1 ? 0 : t+1;
    
    int IarrayNum[3] = {iPrev, i, iNext};
    int JarrayNum[3] = {jPrev, j, jNext};
    int TarrayNum[3] = {tPrev, t, tNext};
    
    
    
    
    for(int sindex=0;sindex<3;sindex++){
        for(int kindex=0; kindex<3;kindex++){
            for(int windex=0; windex<3;windex++){
                
                int s = IarrayNum[sindex];
                int k = JarrayNum[kindex];
                int w = TarrayNum[windex];
                
                int sPrev = s == 0 ? Lx-1 : s-1;
                int sNext = s == Lx-1 ? 0 : s+1;
                int kPrev = k == 0 ? Ly-1 : k-1;
                int kNext = k == Ly-1 ? 0 : k+1;
                int wPrev = w == 0 ? Lt-1 : w-1;
                int wNext = w == Lt-1 ? 0 : w+1;
                
                if(s==i && k==j && w==t){
                    H4a[s][k][w]=H4aflip[s][k][w];
                    H4b[s][k][w]=H4bflip[s][k][w];
                    H4c[s][k][w]=H4cflip[s][k][w];}
                
                else if(sPrev==i && k==j && w==t){
                    H4a[sPrev][k][w]=H4aflip[sPrev][k][w];
                    H4b[sPrev][k][w]=H4bflip[sPrev][k][w];
                    H4c[sPrev][k][w]=H4cflip[sPrev][k][w];}
                
                else if(s==i && kPrev==j && w==t){
                    H4a[s][kPrev][w]=H4aflip[s][kPrev][w];
                    H4b[s][kPrev][w]=H4bflip[s][kPrev][w];
                    H4c[s][kPrev][w]=H4cflip[s][kPrev][w];}
                
                else if(sNext==i && k==j && w==t){
                    H4a[sNext][k][w]=H4aflip[sNext][k][w];
                    H4b[sNext][k][w]=H4bflip[sNext][k][w];
                    H4c[sNext][k][w]=H4cflip[sNext][k][w];}
                
                else if(s==i && kNext==j && w==t){
                    H4a[s][kNext][w]=H4aflip[s][kNext][w];
                    H4b[s][kNext][w]=H4bflip[s][kNext][w];
                    H4c[s][kNext][w]=H4cflip[s][kNext][w];}
                
                
                else if(s==i && k==j && wPrev==t){
                    H4a[s][k][wPrev]=H4aflip[s][k][wPrev];
                    H4b[s][k][wPrev]=H4bflip[s][k][wPrev];
                    H4c[s][k][wPrev]=H4cflip[s][k][wPrev];}
                
                else if(s==i && k==j && wNext==t){
                    H4a[s][k][wNext]=H4aflip[s][k][wNext];
                    H4b[s][k][wNext]=H4bflip[s][k][wNext];
                    H4c[s][k][wNext]=H4cflip[s][k][wNext];}
                
                
                
            
            
            eflip1+= EConfiguration(s,k,w);
            }}}
    
    
    
    return eflip1;}









// -------- Metropolis Update -------
bool UpdateHiggsFieldsRandom () {
    
    int i = int(Lx*std_rand());             // choose a random spin
    int j = int(Ly*std_rand());
	int t = int(Lt*std_rand());
    
    
    //Update H1 amplitude
    
    std::complex<double> Einitial =  E0(i, j, t);
    
    std::complex<double> Efinal   =  EfH1(i, j, t);
    
    if( std_rand() < exp(-1./T*abs(Efinal-Einitial))){
        
        
        H1a[i][j][t]=H1aflip[i][j][t];
        H1b[i][j][t]=H1bflip[i][j][t];
        H1c[i][j][t]=H1cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H1a[i][j][t]=H1a0[i][j][t];
        H1b[i][j][t]=H1b0[i][j][t];
        H1c[i][j][t]=H1c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    //Rotate Rab.H1
 /*
    
    Efinal   =  EfH1Rab(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H1a[i][j][t]=H1aflip[i][j][t];
        H1b[i][j][t]=H1bflip[i][j][t];
        H1c[i][j][t]=H1cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H1a[i][j][t]=H1a0[i][j][t];
        H1b[i][j][t]=H1b0[i][j][t];
        H1c[i][j][t]=H1c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    //Rotate Rbc.H1
    
    
    Efinal   =  EfH1Rbc(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H1a[i][j][t]=H1aflip[i][j][t];
        H1b[i][j][t]=H1bflip[i][j][t];
        H1c[i][j][t]=H1cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H1a[i][j][t]=H1a0[i][j][t];
        H1b[i][j][t]=H1b0[i][j][t];
        H1c[i][j][t]=H1c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    */
    
    //Amplitude H2
    
    
    Efinal   =  EfH2(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H2a[i][j][t]=H2aflip[i][j][t];
        H2b[i][j][t]=H2bflip[i][j][t];
        H2c[i][j][t]=H2cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H2a[i][j][t]=H2a0[i][j][t];
        H2b[i][j][t]=H2b0[i][j][t];
        H2c[i][j][t]=H2c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    //Rotate Rab.H2
    
    
    Efinal   =  EfH2Rab(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H2a[i][j][t]=H2aflip[i][j][t];
        H2b[i][j][t]=H2bflip[i][j][t];
        H2c[i][j][t]=H2cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H2a[i][j][t]=H2a0[i][j][t];
        H2b[i][j][t]=H2b0[i][j][t];
        H2c[i][j][t]=H2c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    //Rotate Rbc.H2
    
    
    Efinal   =  EfH2Rbc(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H2a[i][j][t]=H2aflip[i][j][t];
        H2b[i][j][t]=H2bflip[i][j][t];
        H2c[i][j][t]=H2cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H2a[i][j][t]=H2a0[i][j][t];
        H2b[i][j][t]=H2b0[i][j][t];
        H2c[i][j][t]=H2c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    
    //Amplitude H3
    
    
    Efinal   =  EfH3(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H3a[i][j][t]=H3aflip[i][j][t];
        H3b[i][j][t]=H3bflip[i][j][t];
        H3c[i][j][t]=H3cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H3a[i][j][t]=H3a0[i][j][t];
        H3b[i][j][t]=H3b0[i][j][t];
        H3c[i][j][t]=H3c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    //Rotate Rab.H3
    
    
    Efinal   =  EfH3Rab(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H3a[i][j][t]=H3aflip[i][j][t];
        H3b[i][j][t]=H3bflip[i][j][t];
        H3c[i][j][t]=H3cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H3a[i][j][t]=H3a0[i][j][t];
        H3b[i][j][t]=H3b0[i][j][t];
        H3c[i][j][t]=H3c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    //Rotate Rbc.H3
    
    
    Efinal   =  EfH3Rbc(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H3a[i][j][t]=H3aflip[i][j][t];
        H3b[i][j][t]=H3bflip[i][j][t];
        H3c[i][j][t]=H3cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H3a[i][j][t]=H3a0[i][j][t];
        H3b[i][j][t]=H3b0[i][j][t];
        H3c[i][j][t]=H3c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    
    //Amplitude H4
    
    
    Efinal   =  EfH4(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H4a[i][j][t]=H4aflip[i][j][t];
        H4b[i][j][t]=H4bflip[i][j][t];
        H4c[i][j][t]=H4cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H4a[i][j][t]=H4a0[i][j][t];
        H4b[i][j][t]=H4b0[i][j][t];
        H4c[i][j][t]=H4c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    //Rotate Rab.H4
    
    
    Efinal   =  EfH4Rab(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H4a[i][j][t]=H4aflip[i][j][t];
        H4b[i][j][t]=H4bflip[i][j][t];
        H4c[i][j][t]=H4cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H4a[i][j][t]=H4a0[i][j][t];
        H4b[i][j][t]=H4b0[i][j][t];
        H4c[i][j][t]=H4c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    //Rotate Rbc.H4
    
    
    Efinal   =  EfH4Rbc(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H4a[i][j][t]=H4aflip[i][j][t];
        H4b[i][j][t]=H4bflip[i][j][t];
        H4c[i][j][t]=H4cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H4a[i][j][t]=H4a0[i][j][t];
        H4b[i][j][t]=H4b0[i][j][t];
        H4c[i][j][t]=H4c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
}


bool UpdateHiggsFieldsSequential (int i, int j, int t) {
    

    
    
    //Update H1 amplitude
    
    std::complex<double> Einitial =  E0(i, j, t);
    
    std::complex<double> Efinal   =  EfH1(i, j, t);
    
    if( std_rand() < exp(-1./T*abs(Efinal-Einitial))){
        
        
        H1a[i][j][t]=H1aflip[i][j][t];
        H1b[i][j][t]=H1bflip[i][j][t];
        H1c[i][j][t]=H1cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H1a[i][j][t]=H1a0[i][j][t];
        H1b[i][j][t]=H1b0[i][j][t];
        H1c[i][j][t]=H1c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
 /*   //Rotate Rab.H1
    
    
    Efinal   =  EfH1Rab(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H1a[i][j][t]=H1aflip[i][j][t];
        H1b[i][j][t]=H1bflip[i][j][t];
        H1c[i][j][t]=H1cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H1a[i][j][t]=H1a0[i][j][t];
        H1b[i][j][t]=H1b0[i][j][t];
        H1c[i][j][t]=H1c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    //Rotate Rbc.H1
    
    
    Efinal   =  EfH1Rbc(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H1a[i][j][t]=H1aflip[i][j][t];
        H1b[i][j][t]=H1bflip[i][j][t];
        H1c[i][j][t]=H1cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H1a[i][j][t]=H1a0[i][j][t];
        H1b[i][j][t]=H1b0[i][j][t];
        H1c[i][j][t]=H1c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    */
    
    
    //Amplitude H2
    
    
    Efinal   =  EfH2(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H2a[i][j][t]=H2aflip[i][j][t];
        H2b[i][j][t]=H2bflip[i][j][t];
        H2c[i][j][t]=H2cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H2a[i][j][t]=H2a0[i][j][t];
        H2b[i][j][t]=H2b0[i][j][t];
        H2c[i][j][t]=H2c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    //Rotate Rab.H2
    
    
    Efinal   =  EfH2Rab(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H2a[i][j][t]=H2aflip[i][j][t];
        H2b[i][j][t]=H2bflip[i][j][t];
        H2c[i][j][t]=H2cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H2a[i][j][t]=H2a0[i][j][t];
        H2b[i][j][t]=H2b0[i][j][t];
        H2c[i][j][t]=H2c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    //Rotate Rbc.H2
    
    
    Efinal   =  EfH2Rbc(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H2a[i][j][t]=H2aflip[i][j][t];
        H2b[i][j][t]=H2bflip[i][j][t];
        H2c[i][j][t]=H2cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H2a[i][j][t]=H2a0[i][j][t];
        H2b[i][j][t]=H2b0[i][j][t];
        H2c[i][j][t]=H2c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    
    //Amplitude H3
    
    
    Efinal   =  EfH3(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H3a[i][j][t]=H3aflip[i][j][t];
        H3b[i][j][t]=H3bflip[i][j][t];
        H3c[i][j][t]=H3cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H3a[i][j][t]=H3a0[i][j][t];
        H3b[i][j][t]=H3b0[i][j][t];
        H3c[i][j][t]=H3c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    //Rotate Rab.H3
    
    
    Efinal   =  EfH3Rab(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H3a[i][j][t]=H3aflip[i][j][t];
        H3b[i][j][t]=H3bflip[i][j][t];
        H3c[i][j][t]=H3cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H3a[i][j][t]=H3a0[i][j][t];
        H3b[i][j][t]=H3b0[i][j][t];
        H3c[i][j][t]=H3c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    //Rotate Rbc.H3
    
    
    Efinal   =  EfH3Rbc(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H3a[i][j][t]=H3aflip[i][j][t];
        H3b[i][j][t]=H3bflip[i][j][t];
        H3c[i][j][t]=H3cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H3a[i][j][t]=H3a0[i][j][t];
        H3b[i][j][t]=H3b0[i][j][t];
        H3c[i][j][t]=H3c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    
    //Amplitude H4
    
    
    Efinal   =  EfH4(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H4a[i][j][t]=H4aflip[i][j][t];
        H4b[i][j][t]=H4bflip[i][j][t];
        H4c[i][j][t]=H4cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H4a[i][j][t]=H4a0[i][j][t];
        H4b[i][j][t]=H4b0[i][j][t];
        H4c[i][j][t]=H4c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    //Rotate Rab.H4
    
    
    Efinal   =  EfH4Rab(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H4a[i][j][t]=H4aflip[i][j][t];
        H4b[i][j][t]=H4bflip[i][j][t];
        H4c[i][j][t]=H4cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H4a[i][j][t]=H4a0[i][j][t];
        H4b[i][j][t]=H4b0[i][j][t];
        H4c[i][j][t]=H4c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
    //Rotate Rbc.H4
    
    
    Efinal   =  EfH4Rbc(i, j, t);
    
    if( std_rand() < exp(-1./T*real(Efinal-Einitial))){
        
        
        H4a[i][j][t]=H4aflip[i][j][t];
        H4b[i][j][t]=H4bflip[i][j][t];
        H4c[i][j][t]=H4cflip[i][j][t];
        
        Einitial = Efinal;
        
    }
    else{
        H4a[i][j][t]=H4a0[i][j][t];
        H4b[i][j][t]=H4b0[i][j][t];
        H4c[i][j][t]=H4c0[i][j][t];
        
        Einitial = Einitial;
        
    }
    
    
    
}




// ------- MC Sweep --------------
void MCSweepHiggsFieldsRandom( ) {
    int acceptsH1 = 0, XIaccepts = 0;
    for (int i = 0; i < Lx; i++){
        for (int j = 0; j < Ly; j++){
            for (int t = 0; t < Lt; t++){
                if (UpdateHiggsFieldsRandom())

                    ++steps;
        
            }}}
    }

void MCSweepHiggsFieldsSequential( ) {
    int acceptsH1 = 0, XIaccepts = 0;
    for (int i = 0; i < Lx; i++){
        for (int j = 0; j < Ly; j++){
            for (int t = 0; t < Lt; t++){
                if (UpdateHiggsFieldsSequential(i,j,t))
                    
                    ++steps;
                
            }}}
}



//------- Magnetization function : Qlm ----------



double Q11av() {  //Modulus of magnetisation
    long double Q11Sum = 0;
    int i, j, t;
    
    for (int i = 0; i < Lx; i++){
        for (int j = 0; j < Ly; j++) {
            for (int t = 0; t < Lt; t++) {
                
                
                Q11Sum += H1a[i][j][t]*H1a[i][j][t] + H1b[i][j][t]*H1b[i][j][t] + H1c[i][j][t]*H1c[i][j][t];

                
                
                
            }}}
    
    return double (Q11Sum/(1.*N));
    
}



int main (int argc, char *argv[]) {
    
    
    /*   struct timeval t;
     assert(0 == gettimeofday(&t, NULL));
     srandom(t.tv_usec);*/
    
    srand (time(NULL));
    
    // Get starting timepoint
    auto start = high_resolution_clock::now();
    
    /*cout << " Two-dimensional Ising Model - Metropolis simulation\n"
    << " ---------------------------------------------------\n"
    << " Enter number of spins L in each spatial direction: ";
    cin >> Lx;
    cout << " Enter Initial Temperature T0: ";
    cin >> T0; */


    //------------------ Choose thermodynamic parameters ---------------------------
    Lx=4;
    Ly=4;
    Lt=4;
    
    N = Lx*Ly*Lt;
    
    



    
    //------------------ Choose field theory parameters ---------------------------

    ms=-5.; //WARNING was -5
    u0=1.; //WARNING was 5
    u1=1.; //Warning -1 for U1; +1 for Z2

    
    Keff = 1.;
    Geff = 3.;
    
    g=1.0; //WARNING was g=5
    
    H0=1.0*sqrt(2.*abs(ms))/sqrt(4.*u0 + u1);
    
    
    //------------------ Choose update parameters ---------------------------
    
    d_Phi=0.3*PI;
    
    d_Amp0 = 0.1*H0;
    d_Amp = 0.1*H0;
    
    
    //------------------ Begin First Simulations ---------------------------

    
    int MCSteps=12500;
    
    
    

    InitializeConfiguration();
    
    
    int thermSteps = int(1 * MCSteps);       // thermalization

    
    
    
    // Cooling
    for(int w = 0; w<thermSteps; w++){
        
        MCSweepHiggsFieldsRandom( ); //Random
        MCSweepHiggsFieldsSequential( ); //Sequential
        
    }
    
    cout << "........................................................" << endl;
    cout << "Cooling complete" <<endl;
    cout << "........................................................" << endl;
    
    
    for (int K = 0; K < MCSteps; K++) {
        
        MCSweepHiggsFieldsRandom( ); //Random
        MCSweepHiggsFieldsSequential( ); //Sequential

        double Q11MCav;
        
        Q11MCav += Q11av()/double(MCSteps);
        
        
        
    }
    
    
    //-----------------------------------------------

    // Get ending timepoint
    auto stop = high_resolution_clock::now();
    
    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);
        
    cout << "........................................................" << endl;
    cout << "AcceptanceRatioAx= " << acceptsAx/double(2.*3.*N*MCSteps) <<endl;
    cout << "AcceptanceRatioAy= " << acceptsAy/double(2.*3.*N*MCSteps) <<endl;
    cout << "AcceptanceRatioAxAmp= " << acceptsAxAmp/double(N*MCSteps) <<endl;
    cout << "AcceptanceRatioAyAmp= " << acceptsAyAmp/double(N*MCSteps) <<endl;
    cout << "Time taken to minimize: "<<  duration.count()/pow(10,6)/60 << " minutes" << endl;
    cout << "........................................................" << endl;

}