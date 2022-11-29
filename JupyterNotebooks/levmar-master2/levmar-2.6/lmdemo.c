/////////////////////////////////////////////////////////////////////////////////
// 
//  Demonstration driver program for the Levenberg - Marquardt minimization
//  algorithm
//  Copyright (C) 2004-05  Manolis Lourakis (lourakis at ics forth gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////

/******************************************************************************** 
 * Levenberg-Marquardt minimization demo driver. Only the double precision versions
 * are tested here. See the Meyer case for an example of verifying the Jacobian 
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "levmar.h"





void gauss_f(double* p, double* x, int m, int n, void* data)
{
    double ndbl = n;
    double framesize = sqrt(ndbl);
    //gets data
    int xsize = framesize; //n should be a square number so this is ok
    int ysize = xsize;  //assumes square data area


    //grabs the current guesses
    int p_count = 0;
    double A = p[p_count++];
    double sigma_xy = p[p_count++];
    double b = p[p_count++];
    double Xo = p[p_count++];
    double Yo = p[p_count++];

    //initialises other useful variables
    double X, Y, e;
    int element;

    //loops over all pixels in the image
    for (int y_coord = 0; y_coord < ysize; y_coord++)
    {
        for (int x_coord = 0; x_coord < xsize; x_coord++)
        {
            X = x_coord + 1; //+1 because we start at 0 in the for loop
            //but the real data will be from a matrix indexed to start at [1,1]
            Y = y_coord + 1;

            element = x_coord + (y_coord * ysize);
            //y_coord *y_size so it starts at 0, then 1*ysize...so element[2,1] in the original matrix is
            //element y_size+1 in f etc

            //=======Model Function - amplitude * 2D gaussian + background==========
            e = exp(-((pow(X - Xo, 2) + pow(Y - Yo, 2)) * 0.5 * pow(sigma_xy, -2))); //2D gaussian
            x[element] = (A * e) + b;
            //======================================================================
        }
    }



}

void gauss_jac(double* p, double* jac, int m, int n, void* data)
//calculates the Jacobian vector of derivatives for each of the variables to be fitted

{
    double ndbl = n;
    double framesize = sqrt(ndbl);
    //gets data
    int xsize = framesize; //n should be a square number so this is ok
    int ysize = xsize;  //assumes square data area




    //grabs the initial guesses
    int p_count = 0;
    double A = p[p_count++];
    double sigma_xy = p[p_count++];
    double b = p[p_count++];
    double Xo = p[p_count++];
    double Yo = p[p_count++];
    //initialises other useful variables
    double X, Y, e;
    int element;
    int j = 0;

    //loops over all pixels in the image
    for (int y_coord = 0; y_coord < ysize; y_coord++)
    {
        for (int x_coord = 0; x_coord < xsize; x_coord++)
        {
            // Jacobian matrix J(i,j) = dfi / dxj,
            // fi = 2d gaussian with background
            // and the xj are the parameters (A,sigma,b,Xo,Yo)

            X = x_coord + 1;
            Y = y_coord + 1;
            element = x_coord + (y_coord * ysize);
            e = exp(-((pow(X - Xo, 2) + pow(Y - Yo, 2)) * 0.5 * pow(sigma_xy, -2)));



            jac[j++] = e; //dF/dA
            ///can put if variant =! 2 then run this next line style thing
            jac[j++] = pow(sigma_xy, -3) * (pow(X - Xo, 2) + pow(Y - Yo, 2)) * A * e;  //dF/d(sigma)
            jac[j++] = 1.0;                                    //dF/db
            jac[j++] = pow(sigma_xy, -2) * (X - Xo) * A * e; //dF/dXo
            jac[j++] = pow(sigma_xy, -2) * (Y - Yo) * A * e; //dF/dYo
        }
    }

}


int main()
{

    double p[5] = { 5, 1, 5, 2, 2 };
    double x[16] = { 10, 10, 15, 10,
            10, 15, 20, 15,
            10, 10, 15, 10,
            10, 10, 10, 10 };
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0] = LM_INIT_MU;
    opts[1] = 0.001;
    opts[2] = 0.00001;
    opts[3] = 0.001;
    opts[4] = LM_DIFF_DELTA;
    int m, n;
    m = 5;
    n = 16;
    int ret;
    //double lb[5] = {0,0,0,0,0};
    //double ub[5] = { DBL_MAX, 5, DBL_MAX, 5, 5 };
    //double dscl[3];
    ret = dlevmar_der(gauss_f,  //pointer to the function that calculates the model values at each datapoint given parameters p
            gauss_jac, //pointer to the function that calculates the analytic jacobian dxi/dpj at each datapoint
            p, //current best guess parameters as an array of size m
            x,  //data to be fitted to the model, array of size n
            m,  //number of free parameters
            n, //number of datapoints

            ///constraints
            //lb, //array of m elements setting lower bounds on each parameter (p[i] >= lb[i])
            //ub, //array of m elements setting upper bounds on each parameter (p[i] <= ub[i])
            //dscl, //???
            1000, //max no. of iterations the solver will run for, 50 to 100 seems about right
            opts, //5 element array, described above
            info, //from source code of levmar library:
            /* info: information regarding the minimization. Used to find the lsq minimisation value at output
      * info[0]= ||e||_2 at initial p.
      * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
      * info[5]= # iterations,
      * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
      *                                 2 - stopped by small Dp
      *                                 3 - stopped by itmax
      *                                 4 - singular matrix. Restart from current p with increased mu
      *                                 5 - no further error reduction is possible. Restart with increased mu
      *                                 6 - stopped by small ||e||_2
      *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
      * info[7]= # function evaluations
      * info[8]= # Jacobian evaluations
      * info[9]= # linear systems solved, i.e. # attempts for reducing error
      */
      
            NULL,    // working memory at least LM_BLEC_DER_WORKSZ() reals large, allocated if NULL
            NULL,   //Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed.
            NULL); // pointer to additional data passed to function and jacobian calculator - fixed constants etc
            


    printf("Results for Gaussian fitting:\n");
    printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
    register int i;
    for (i = 0; i < m; ++i)
        printf("%.7g ", p[i]);
    printf("\n\nMinimization info:\n");
    for (i = 0; i < LM_INFO_SZ; ++i)
        printf("%g ", info[i]);
    printf("\n");

  
    printf("New version\n");

    return 0;
}
