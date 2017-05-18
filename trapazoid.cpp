//===Project: DualPol ===== Author: John M. Krause =====Nov. 2001 =============
//
// Module Name: trapazoid.cpp 
//
// Module Version: 1.0
//
// Module Language: c++
//
// Change History:
//
//Date    Version    Programmer  Notes
//-------------------------------------------------------------
//06/19/02  1.0    John Krause  Creation
//
//
// Calling Sequence: 
// 
//
// Description:
// This subroutine calcuates the value of the defined trapazoid 
//
//
//(c) Copyright 2004, Board of Regents of the University of Oklahoma.
//All rights reserved. Not to be provided or used in any format without
//the express written permission of the University of Oklahoma.
//Software provided as is no warranties expressed or implied.
//========================================================================
//Specification Files:
#include <cassert>
#include <cstdlib> // exit()
#include <cmath> // finite()
#include <iostream>
#include <dualpol/trapazoid.h>

using std::cout;
using std::endl;

//========================================================================

float trap( float input, float low, float high, float interval)
{
  if ( !std::isfinite(input)) {
      cout << "Problem with input value in trapazoid calculation: " << endl
      << " Value: " << input << endl;
       return 0.0;
  }

  if ( input >= low && input <= high )
    return 1.0;
  if ( input < (low-interval))
    return 0.0;
  if ( input > (high+interval))
    return 0.0;
  if ( input >= (low-interval) && input < low )
    return (1.0 - (low - input)/interval);
  if ( input <= (high+interval) && input > high )
    return (1.0 - (input-high)/interval);

  cout << "Problem with trapazoid calculation: " << endl
       << " Low: " << low << endl
       << " High: " << high << endl
       << " Interval: " << interval << endl
       << " Value: " << input << endl;
  exit(-1);
}

//====end of trapazoid ====================================================
float trap4point( float input, float x1, float x2, float x3, float x4)
{
  //abort for poor membership functions
  if ( (x2-x1) < 0 || (x3-x2) < 0 || (x4-x3) < 0 ) return 0;
  //abort for poor input values
  if ( !std::isfinite(input)) {
      cout << "Problem with input value in trapazoid calculation: " << endl
      << " Value: " << input << endl;
       return 0.0;
  }

  assert(x1<=x2);
  assert(x2<=x3);
  assert(x3<=x4);

  if ( input >= x2 && input <= x3 )
    return 1.0;
  if ( input <= x1 )
    return 0.0;
  if ( input >= x4 )
    return 0.0;
  if ( input > x1 && input < x2 )
    return ( (input-x1)/ (x2-x1));
  if ( input > x3 && input < x4 )
    return ( (x4-input)/(x4-x3) );

  cout << "Problem with trapazoid calculation:" << endl
       << " x1: " << x1 << endl
       << " x2: " << x2 << endl
       << " x3: " << x3 << endl
       << " x4: " << x4 << endl
       << " Value: " << input << endl;
  exit(-1);
}
//====end of trapazoid ====================================================
