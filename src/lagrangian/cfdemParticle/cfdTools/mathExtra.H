/* ----------------------------------------------------------------------
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
    
    Copyright of this contribution:
    Copyright 2014-     TU Graz, IPPT
------------------------------------------------------------------------- */

#ifndef CFDEM_MATH_EXTRA_H
#define CFDEM_MATH_EXTRA_H

#include <vector>
//#include "math.h"
#include "stdio.h"
#include "string.h"
#include "error.h"
#include "ctype.h"

#define TOLERANCE_ORTHO 1e-10

namespace MathExtra 
{

//  inline void     outerProduct(double *vec1, double *vec2, double **m);
//  inline double  spheroidGeometry(int index, double bi, double ai);



//--------------------------------------------------------------------
//   Outer Product of two vectors
inline void outerProduct(double *vec1, double *vec2, double **m)
{
  int i, j;
  //debug output
//  for( i = 0; i < 3; ++i )
//    printf("OUTER PRODUCT: Input: vec1 element %d = %g", i, vec1[i]);
//  for( i = 0; i < 3; ++i )
//    printf("OUTER PRODUCT: Input: vec2 element %d=%g", i, vec2[i]);
  
  //calculation
  for( i = 0; i < 3; ++i )
    for( j = 0; j < 3; ++j )
    {
      m[i][j] = vec1[i] * vec2[j];
      printf("OUTER PRODUCT: Result: m[%d][%d]=%g", i, j, m[i][j]);
    }
    
    
}

//--------------------------------------------------------------------
// Compute the diameter of the hydrodynamically equivalent cylinder
inline double spheroidDiameter(
                                double shapeX, double shapeY, double shapeZ,
                                double& aspectRatio
                              )
{
    //INPUT
    //  shape       ...values with the half-axes of the spheroid
    
    //OUTPUT
    //  aspectRatio ...the aspect ratio of the spheroid (ai/bi)
    //                 is also that of the cylinder element
    //  return value... the hydrodynamically equivalent diameter 

     double   bi = std::min( shapeX, std::min(shapeY,shapeZ) );
     aspectRatio = std::max( shapeX, std::max(shapeY,shapeZ) ) 
                 / std::max(bi,1e-32);
     aspectRatio = std::max(1.000001,aspectRatio);

     return 2.48 * bi / std::sqrt( std::log(aspectRatio) );

}

//--------------------------------------------------------------------
// Compute the major, minor axis and eccentricity parameters of a prolate spheroid
inline bool spheroidGeometry(double aspectRatio,     //input
                             double& ei, double& Le              //outputs
                            ) //
{
    //INPUT
    //  aspectRatio ...major/minor aspect ratio
    
    //OUTPUT
    //  ei  ... 
    //  Le  ... 

 	ei = std::sqrt( 
                     1.0
                   - 1.0 / (aspectRatio*aspectRatio)
 	              );
 	Le = std::log(
 		           (1.0+ei)
 		          /(1.0-ei)
                 );

	return true; 		    
}

//--------------------------------------------------------------------
// Compute the major, minor axis and eccentricity parameters of a prolate spheroid
inline double Pi()
{
    return 3.1415926535897932384626433832795;
}


//--------------------------------------------------------------------
// Compute the eccentricity parameters of a prolate spheroid
inline bool spheroidGeometry2(   double  aspectRatio,                 //inputs
                                 double& XAe,    double& YAe,         //outputs
                                 double& XCe,    double& YCe,         //outputs
                                 double& YHe                          //outputs
                                ) //
{
    //INPUT
    //  aspectRatio ...major/minor aspect ratio of the spheroid
    
    //OUTPUT
    //  XAe  ...Eccentricity dependet parameter
    //  YAe  ...Eccentricity dependet parameter
    //  XCe  ...Eccentricity dependet parameter
    //  XCe  ...Eccentricity dependet parameter
    //  YCe  ...Eccentricity dependet parameter
    //  YHe  ...Eccentricity dependet parameter
  

    double ei(0.0), Le(0.0);  
    bool   result = 
           spheroidGeometry(aspectRatio,  //inputs
                            ei,    Le     //outputs
                           );
    if(!result)
        return false;
        
    XAe= 2.6666666666666666666666667
         *ei*ei*ei
         /(-2.0*ei+(1.0+ei*ei)*Le); 
    YAe= 5.333333333333333333333333333
        *ei*ei*ei
        /(2.0*ei+(3*ei*ei-1.0)*Le);
    XCe= 1.333333333333333333333333333
        *ei*ei*ei
        *(1.0-ei*ei)
        /(2.0*ei-(1.0-ei*ei)*Le);
    YCe= 1.3333333333333333333333
        *ei*ei*ei
        *(2.0-ei*ei)
        /(-2.0*ei+(1.0+ei*ei)*Le);
    YHe= 1.3333333333333333333333
        *ei*ei*ei*ei*ei
        /(-2.0*ei+(1.0+ei*ei)*Le);

	return true; 		    

}

//--------------------------------------------------------------------
// zeroize a 3x3x3 tensor
inline void zeroize333(double tensor[3][3][3] )
{
    for(int iX=0; iX<3; iX++)
        for(int iY=0; iY<3; iY++)
            for(int iZ=0; iZ<3; iZ++)
                tensor[iX][iY][iZ] = 0.0;
}

//--------------------------------------------------------------------
// zeroize a 3x3 tensor
inline void zeroize33(double tensor[3][3] )
{
    for(int iX=0; iX<3; iX++)
        for(int iY=0; iY<3; iY++)
                tensor[iX][iY] = 0.0;
}

//--------------------------------------------------------------------
// multiply a 3x3x3 tensor with a scalar
inline void multiply333(double scalar, double tensor[3][3][3] )
{
    for(int iX=0; iX<3; iX++)
        for(int iY=0; iY<3; iY++)
            for(int iZ=0; iZ<3; iZ++)
                tensor[iX][iY][iZ] *= scalar;
}
//--------------------------------------------------------------------
// Compute a dot and dyadic product of with a vector
inline void permutationTensor(double tensor[3][3][3] )
{
    zeroize333(tensor);
    tensor[0][1][2] = 1.0; 
    tensor[1][2][0] = 1.0; 
    tensor[2][0][1] = 1.0; 
    tensor[0][2][1] =-1.0; 
    tensor[2][1][0] =-1.0; 
    tensor[1][0][2] =-1.0; 
}



//--------------------------------------------------------------------
// Compute a dot product of the permutation tensor and 
// then a dyadic product of with a vector
inline bool permutationDotDyadic( 
                                 double vector[3], 
                                 double  tensor[3][3][3]      
                                )
{
    //Generate permutation tensor
    double permutationT[3][3][3];
    permutationTensor(permutationT);
    
    //Step 1: compute dot prodcut of permutation tensor
    double permutationDotProd[3][3];
    zeroize33(permutationDotProd);

    for(int iX=0; iX<3; iX++)
        for(int iY=0; iY<3; iY++)
            for(int iZ=0; iZ<3; iZ++)
                permutationDotProd[iX][iY] += permutationT[iX][iY][iZ]
                                            * vector[iZ];

    for(int iX=0; iX<3; iX++)
        for(int iY=0; iY<3; iY++)
            for(int iZ=0; iZ<3; iZ++)
                tensor[iX][iY][iZ] = permutationDotProd[iX][iY]
                                   * vector[iZ]; 
    return true;

}

//--------------------------------------------------------------------
// Compute a dot and dyadic product of with a vector
inline bool doubleDotTensor333Tensor33(double tensor333[3][3][3], 
                                       double tensor33[3][3],
                                       double result[3]
                                      )
{
    result[0]=0.0;result[1]=0.0;result[2]=0.0;
    
    for(int iX=0; iX<3; iX++)
        for(int iY=0; iY<3; iY++)
            for(int iZ=0; iZ<3; iZ++)
                result[iX] += tensor333[iX][iY][iZ] * tensor33[iY][iZ];

    return true;
}

}; //end of namespace

#endif
