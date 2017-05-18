//===Project: DualPol ===== Author: John M. Krause ===== April 2002 ===========
//
// Module Name: calc_texture.cpp 
//
// Module Version: 1.0
//
// Module Language: c++
//
// Change History:
//
//Date    Version    Programmer  Notes
//-------------------------------------------------------------
//11/08/02  1.0    John Krause  Creation, as calc_STD_Z
//04/09/04  1.1    John Krause  Generic texture calcs now 
//12/10/04	1.2 	 John Krause	Calc_texture now uses biased STD from 
//                              standard deviation
//
//
// Calling Sequence: 
// 
//
// Description:
//
//
//
//(c) Copyright 2004, Board of Regents of the University of Oklahoma.
//All rights reserved. Not to be provided or used in any format without
//the express written permission of the University of Oklahoma.
//Software provided as is no warranties expressed or implied.
//========================================================================
//Specification Files:
#include <cmath> // fmod
#include <cstring> //memcpy 
#include <iostream>
#include <vector>
#include <dualpol/Average_Filter.h>
#include <dualpol/DualPolConsts.h> // NUM_GATES
#include <dualpol/calc_texture.h>
#include <dualpol/standard_deviation.h>
//========================================================================
void calc_texture( const int num_bins,
                   const float no_data_value,
                   const float input[],
                   float STD_output[],
                   const float std_thresholds,
                   const int filter_length )
{
  int i;
  float smooth_input[num_bins];
  float no_data = no_data_value;

//first the Filter length should be an odd number so that we can appropriately
//center the value.

  if ( fmod(filter_length, 2.0 ) == 0 ) {
    std::cout << "Warning!! Filtering for an even Filter Length is not "
              << "recommended. Please select an odd Filter Length "
              << filter_length << std::endl;
  }
//
  int center = (int) floor( filter_length/2.0 );
  int start = center;
  int end = num_bins-center; 
  std::vector<float> diff_data;
//
// copy Z 
//
/*
  for(i=0;i<num_bins;i++) {
    smooth_input[i] = input[i]; 
  }
*/
//Hoyt suggest this to improve speed of execution 
//12-18-2014 --jmk
   memcpy(smooth_input, input, (sizeof(float)*num_bins));
//
// apply an Average Filter to input to get smoothed input
//
  if (filter_length >= 0 ) {   
    Average_Filter( filter_length, num_bins, no_data, smooth_input);
  } else {
    std::cout << " no average filtering applied to calc_texture: \n";
    std::cout << " you may have issues with this...: " << filter_length << "\n";
  }
 
  for(i=0;i<num_bins;i++) {
//
// Calculate the residuals for the standard deviation calculation.
// The residuals tell us how "choppy" Z is and the standard deviation of
// them tells us how extreme that "choppyness" is. If the standard deviation 
// is high then the Z data is bouncing around alot. If low then just a little.
// The point of all this is that STD_Z of AP >> STD_Z of Birds/Insects 
// >> STD_Z of meteorological targets 
//

      start = i-center;
      end = i+center;
   
      if(start < 0 ) start = 0;
      if(end > num_bins-1) end = num_bins-1;

      diff_data.clear();
      for (int j=start;j<=end;j++) {
        if ( input[j] != no_data && 
             smooth_input[j] != no_data && 
             (int) diff_data.size() != filter_length) {

          if ( (input[j]-smooth_input[j]) > std_thresholds || (input[j]-smooth_input[j]) < (-1.0*std_thresholds) || !finite(input[j]) || !finite(smooth_input[j]) ) {
          //this is to sanity check the data. 
          //the data is not always as reliable as we would like it to be
          //you may wish the code to exit if it fails this sanity checking.
          //we choose to plow forward
          //-jmk 12/10/2004
          /*
            cout << "texture: data problemo: " << j;
            cout << " input: " << input[j] << " smo: " << smooth_input[j]; 
            cout << " size: " << diff_data.size() << " filter length: ";
            cout << filter_length << "\n";
          */
          } else {
            diff_data.push_back(input[j]-smooth_input[j]);
          }
        }
      }
//
//check the size of the data
//    
      if ( (int) diff_data.size() >= (int) floor(filter_length/2.0) ) { 
        STD_output[i] = standard_deviation( &diff_data.front(), diff_data.size() );
/*
        cout << "input i: " << i; 
        for (int j=start;j<=end;j++) {
          cout << " " << input[j];
        } 
        cout << "\n";
        cout << "smoothed i: " << i; 
        for (int j=start;j<=end;j++) {
          cout << " " << smooth_input[j];
        } 
        cout << "\n";

        cout << "i: " << i << " STD: " << STD_output[i] << endl;
        for(int k=0;k<(int) diff_data.size();k++)
           cout << "  " << diff_data[k];
        cout << endl;
*/

      } else {
        //cout << "calc_STD_Z: Filter length requirement not meet. \n";
        STD_output[i] = no_data;
      }

    if ( !finite(STD_output[i]) ) {
      std::cout << "There is a problem with the texture calc: " << i << "\n"
                << " STD_Z[i]: " << STD_output[i]
                << " diff_data: \n";
      for(int j=0;j<(int) diff_data.size();j++) {
        std::cout << " j: " << j << "  " << input[i] << " " << smooth_input[j]
                  << "  " << diff_data[j] << "\n";
      }
      std::cout << "Trying to recover.....\n";
      STD_output[i] = no_data;
    }
  }//for i < num_bins 

}
//====end of calc_texture ====================================================

