//===Project: DualPol ===== Author: John M. Krause =====Sept 2012 =============
//
// Module Name: findMetSignal.cpp 
//
// Module Version: 1.0
//
// Module Language: c++
//
// Change History:
//
//Date    Version    Programmer  Notes
//-------------------------------------------------------------
//09/04/12  1.0    John Krause  Creation
//
// Calling Sequence: 
// 
//
// Description:
// Algorithm attempts to find the meteorlogical signal from the data presented.
//
//
//(c) Copyright 2012, Board of Regents of the University of Oklahoma.
//All rights reserved. Not to be provided or used in any format without
//the express written permission of the University of Oklahoma.
//Software provided as is no warranties expressed or implied.
//========================================================================
//Specification Files:
//
#include <algorithm>  //sort
#include <cassert>
#include <cfloat> // FLT_MAX
#include <cmath> // fabs()
#include <cstring> // memcpy()
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

#include <dualpol/trapazoid.h> 
#include <dualpol/calc_texture.h> 
#include <dualpol/utils.h> 
#include <dualpol/Median_Filter.h> 
#include <dualpol/Average_Filter.h> 

//Infrastructure
//
//

using std::cout;
using std::endl;
//========================================================================
void findMetSignal(const int num_bins, 
                     const float no_data, 
                     const float *in_Z, 
                     const float *in_V, 
                     const float *in_SPW, 
                     const float *in_RhoHV, 
                     const float *in_Zdr, 
                     const float *in_PhiDP,
                     const float *in_SNR,
                     float* MetSig) {

int i;
float signal_value = 0.0;
float signal_strength = 0.0;
int num_s = 0;

float *std_PhiDP = new float[num_bins];
float *std_Zdr = new float[num_bins];
float *std_Rho = new float[num_bins];
float *Z = new float[num_bins];
float *V = new float[num_bins];
float *Zdr = new float[num_bins];
float *RhoHV = new float[num_bins];
float *PhiDP = new float[num_bins];
float *SPW = new float[num_bins];
float *SNR = new float[num_bins];


 memcpy(Z, in_Z, sizeof(float)*num_bins);
 memcpy(V, in_V, sizeof(float)*num_bins);
 memcpy(SPW, in_SPW, sizeof(float)*num_bins);
 memcpy(RhoHV, in_RhoHV, sizeof(float)*num_bins);
 memcpy(Zdr, in_Zdr, sizeof(float)*num_bins);
 memcpy(PhiDP, in_PhiDP, sizeof(float)*num_bins);
 memcpy(SNR, in_SNR, sizeof(float)*num_bins);

float weight = 0;

float z_weight = 2.0;
float rhv_weight = 1.0;
float v_weight = 1.0;
float std_zdr_weight = 2.0;
float std_phi_weight = 2.0;
float std_rho_weight = 1.0;

bool output_flag = false;
//compute some additional data
calc_texture(num_bins, no_data, PhiDP, std_PhiDP, 100.0, 9);
calc_texture(num_bins, no_data, Zdr, std_Zdr, 100.0, 9);
calc_texture(num_bins, no_data, RhoHV, std_Rho, 100.0, 9);

for(i=0;i<num_bins;i++) {
  //for safety:
  MetSig[i] = 0.0;

  signal_value = 0.0; //final value
  signal_strength = 0.0; //signal strength
  weight = 0; //weight of those signals
  num_s = 0; //raw number of signals present

  if (output_flag) {
  std::cout << "i: " << i;
  }

  //note that the signal value of Z is doubled
  if ( Z[i] != no_data && std::isfinite(Z[i]) ) {
    signal_value += z_weight*trap4point(Z[i],10.0,30.0,FLT_MAX,FLT_MAX);
    num_s++;
    weight += z_weight;

    if (output_flag) {
      std::cout << " Z: " << Z[i] << " " << trap4point(Z[i],10.0,30.0,FLT_MAX,FLT_MAX);
    }
  }

  //
  //
  if ( RhoHV[i] != no_data && std::isfinite(RhoHV[i]) ) {
    signal_value += rhv_weight*trap4point(RhoHV[i],0.75,0.9,FLT_MAX,FLT_MAX);
    num_s++;
    weight += rhv_weight;

    if (output_flag) {
      std::cout << " Rho: " << RhoHV[i] << " " <<  trap4point(RhoHV[i],0.75,0.9,FLT_MAX,FLT_MAX); 
    }
  }

//
  if ( V[i] != no_data && std::isfinite(V[i]) ) {
  //velocity is an inverted trapazoid. Small values of velocity mean the data
  //are not meteorlogical, while large values mean that the data are.
    signal_value += 1.0 - v_weight*trap4point(V[i],-1.5,-1.0,1.0,1.5);
    num_s++;
    weight += v_weight;

    if (output_flag) {
      std::cout << " V: " << V[i] << " " << (1.0 - trap4point(V[i],-1.5,-1.0,1.0,1.5)); 
    }
  }
//
//We can also derive some parameters from the data which will help
//
//Standard Deviation of PhiDP. Raw values of PhiDP are not useful, but
//PhiDP becomes much "smoother" in the meteorlogical data and the measure of
//this is that the standard deviation of PhiDP is small.
//
//
//compute the STD_PhiDP
  if ( std_PhiDP[i] != no_data && std::isfinite(std_PhiDP[i]) ) {
//a set of values closer to HydroClass.cpp
//    signal_value += std_phi_weight*trap4point(std_PhiDP[i],0.0,1.0,15.0,30.0);
    signal_value += std_phi_weight*trap4point(std_PhiDP[i],0.0,0.0,10.0,20.0);
    num_s++;
    weight += std_phi_weight;

    if (output_flag) {
      std::cout << " std_PhiDP: " << std_PhiDP[i] << " " << trap4point(std_PhiDP[i],0.0,0.0,10.0,20.0); 
    }
  }

//compute the STD_ZDr.
  if ( std_Zdr[i] != no_data && std::isfinite(std_Zdr[i]) ) {
    signal_value += std_zdr_weight*trap4point(std_Zdr[i],0.0,0.0,1.0,2.0);
//a set of values closer to HydroClass.cpp
//    signal_value += std_zdr_weight*trap4point(std_Zdr[i],0.0,0.5,3.0,6.0);
    num_s++;
    weight += std_zdr_weight;
 
    if (output_flag) {
      std::cout << " std_Zdr: " << std_Zdr[i] << " " << trap4point(std_Zdr[i],0.0,0.0,1.0,2.0); 
    }
  }
//compute the STD_Rho_
  if ( std_Rho[i] != no_data && std::isfinite(std_Rho[i]) ) {
    signal_value += std_rho_weight*trap4point(std_Rho[i],0.0,0.0,0.02,0.04);
    num_s++;
    weight += std_rho_weight;

    if (output_flag) {
      std::cout << " std_Rho: " << std_Rho[i] << " " << trap4point(std_Rho[i],0.0,0.0,0.02,0.04); 
    }
  }
//
//
//Locations with large amounts of Texture are non-meteorlogical while areas 
//with moderate to low amounts of texture are meteorlogical
//
//

  //prevents divide by zero errors
  if (num_s >= 4 ) {
    signal_strength = signal_value/weight;

    //scale the Met Signal into a precentage
    MetSig[i] = nint(signal_strength*100.0);

  } else {
      MetSig[i] = -10.0; // = -1.0: meaning I Don't Know? or no value?
  }

  //
//Hardlimits
  //add some hard thresholds:
  //SNR threshold

float  threshold = 80.0;

  if ( SNR[i] < 0.0 ) {
//    MetSig[i] = -10.5;
  }
// RhoHV < 0.65
  if (RhoHV[i] < 0.65 && MetSig[i] >= threshold ) {
    MetSig[i] = -5.0;
  }
//large Zdr
  if (Zdr[i] < -4.5 && MetSig[i] >= threshold ) {
    MetSig[i] = -5.1;
  }
  if (Zdr[i] > 4.5 && MetSig[i] >= threshold ) {
    MetSig[i] = -5.2;
  }

  if (output_flag) {
    cout << " MetSig: " << MetSig[i] << "\n";
  }

/*
  cout << "    i: " << i << " ms: " << MetSig[i] << " z: " << Z[i]; 
  cout << " r: " << RhoHV[i] << " v: " << V[i] << " sp: " << std_PhiDP[i]; 
  cout << " szd: " << std_Zdr[i] << " sr: " << std_Rho[i] << "\n"; 
 */  

  if (MetSig[i] < -10.5 ) {
    cout << "    i: " << i << " ms: " << MetSig[i] << " z: " << Z[i]; 
  }
}



//cleanup
delete [] std_PhiDP;
delete [] std_Zdr;
delete [] std_Rho;

delete [] Z;
delete [] V;
delete [] Zdr;
delete [] RhoHV;
delete [] PhiDP;
delete [] SPW;
delete [] SNR;
}
//========================================================================
