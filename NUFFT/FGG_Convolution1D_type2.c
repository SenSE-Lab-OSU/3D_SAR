/*You can include any C libraries that you normally use*/
#include "math.h"
#include "mex.h"   /*--This C library is required*/

/*
This code implements the expensive convolution loop (see page 448) in
[1] L. Greengard and J.-Y. Lee, "Accelerating the Nonuniform Fast Fourier 
 Transform," SIAM Review, 2004. It is slightly different because we are 
 * doing a type-II transform (uniform --> nonuniform)
 
This code is free to use, but we ask that you please reference the source,
as this will encourage future funding for more free AFRL products. This
code was developed through the AFOSR Lab Task "Moving-Target Radar Feature
Extraction."
Project Manager: Arje Nachman
Principal Investigator: Matthew Ferrara
Date: November 2008
 
Code by (send correspondence to):
Matthew Ferrara, Research Mathematician
AFRL Sensors Directorate Innovative Algorithms Branch (AFRL/RYAT)
Matthew.Ferrara@wpafb.af.mil
 
NOTE: things could be sped up by doing everything with single-precision 
 arithmetic when M_sp<=6 because then the NUFFT is only accurate up to 6 
 digits. However, I have only been able figure out how to pass double 
 precision data to mex functions...
 */

#define PI 3.141592653589793

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*---Inside mexFunction---
    Matlab use:   [f_taur,f_taui]=...
    FGG_Convolution2D_type2(double(real(f_tau(:))),double(imag(f_tau(:))),...
    double(knots),E_3x,[M_sp, tau(1), M_r(1)]);
	*/
    
    /*Declarations*/
    mxArray *ftaupointr;/*The pointer to the real components of the 
                     frequency-domain data*/
    mxArray *ftaupointi;/*The pointer to the imaginary components of the 
                     frequency-domain data*/
    double *ftaur, *ftaui;/*The real and imaginary components of the 
                    frequency-domain data*/
    mxArray *knotspoint;/*The k-space locations*/
    double *knots;
    double *E_2x;
    
    mxArray *E_3xpoint;/*The constant factors of the Gaussian spreading 
                       function*/
    double *E_3x;
    
    mxArray *Scales_point;/*The spreading distance (have far to truncate 
                          the Gaussian)*/
    double *Scales;
    
    double *outArrayR, *outArrayI;
    int i, j, l1, l2;/*Indexing variables*/
    
    double  knotx,  x,  E_1x,  E_2xdummy,  V0r,V0i, 
    V1r, V1i,taux,  M_rx,  M_rxd2, E_23x, E_2xdummy_inv;
    int M, N2, m1, lx, rx, xind, ind, M_sp,  TwoM_sp, R;
    
    /*Copy input pointer fpointr*/
    ftaupointr = prhs[0];
    /*Copy input pointer fpointi*/
    ftaupointi = prhs[1];
    /*Copy input pointer knotspoint*/
    knotspoint = prhs[2];
    //Copy input pointer E_3xspoint
    E_3xpoint = prhs[3];
    /*Copy input pointer Scales_point*/
    Scales_point = prhs[4];
    /*Get vector f*/
    ftaur = mxGetPr(ftaupointr);
    ftaui = mxGetPr(ftaupointi);
    M =  mxGetM(knotspoint);/*number of data points*/
    /*Get matrix knots*/
    knots = mxGetPr(knotspoint);
    /*Get E_3x and E_3y*/
    E_3x = mxGetPr(E_3xpoint);
    /*Get NUFFT parameters*/
    Scales= mxGetPr(Scales_point);
    M_sp=Scales[0];/*Get M_sp*/
    taux=Scales[1];/*Get tau_x*/
    M_rx=Scales[2];/*desired time-domain grid length...assume N is even*/
    /*Store some useful values*/
    N2=M_rx;
    M_rxd2=M_rx/2;
    TwoM_sp=2*M_sp;
    E_2x=(double*)malloc(TwoM_sp*sizeof(double));
    /*Allocate memory and assign output pointer*/
    plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);/*mxReal is our data-type*/
    plhs[1] = mxCreateDoubleMatrix(M, 1, mxREAL);/*mxReal is our data-type*/
    /*Get a pointer to the data space in our newly allocated memory*/
    outArrayR = mxGetPr(plhs[0]);
    outArrayI = mxGetPr(plhs[1]);/*outArrayR+sqrt(-1)*outArrayI = tau =
                                 regularly-sampled Fourier data*/
    /*Initialize the output arrays (M_rx x M_ry matrices)*/
    for(i=0;i<M;i++)
    {
        outArrayR[i]=0;
        outArrayI[i]=0;
    }
    E_2x[M_sp-1]=1;
    /*Perform the (approximate) convolution loop between data and Gaussian*/
    for(i=0;i<M;i++)
    {
        /*store the ith datum's knot location*/
        knotx= knots[i];/*The ith knot's x location*/
        /*Determine the closest index [m1,m2]*/
        m1=floor(M_rx*knotx/(2*PI));
        /*compute the Gaussian factors (we would precompute these in an 
        iterative scheme to save FLOPs) We could compute some of these scalars 
		outside of the loop, too.*/
        x=knotx-m1*PI/M_rxd2;
        E_1x=exp(-x*x/(4*taux));
        E_2xdummy=exp(x*PI/(M_rx*taux));
        E_2xdummy_inv=1/E_2xdummy;
        /*Compute the E_2 vector of powers of the exponential:*/
        /* for(j=0;j<TwoM_sp;j++)//E_2x=E_2xdummy.^((1-M_sp):M_sp);
        // {//Here, we are avoiding unnecessary exponential calculations
        //     E_2x[j] = pow(E_2xdummy,j-M_sp+1);
        }*/
        /*Is pow() the fastest thing we have at our disposal?!? It seems 
        very slow. Instead of the pow() loop, we can calculate:*/
        for(j=M_sp;j<TwoM_sp;j++)/*E_2x=E_2xdummy.^((1-M_sp):M_sp);*/
        {/*Here, we are avoiding unnecessary exponential calculations AND 
         calls to pow()*/
        	E_2x[j] = E_2xdummy*E_2x[j-1];
        }
        for(j=M_sp-2;j>=0;j--)/*E_2x=E_2xdummy.^((1-M_sp):M_sp);*/
        {/*Here, we are avoiding unnecessary exponential calculations AND 
         calls to pow()*/
        	E_2x[j] = E_2x[j+1]*E_2xdummy_inv;
        }
        /*The to small loops above appear to be much faster than the call 
        to pow().  As a non-C programmer, I am very disappointed that I 
        had to do this explicitly.*/
        /*Now we can compute the first constants (see the algorithm 
        description on page 448 of [1])*/
        V0r=E_1x;
        /*Compute f[i]'s contribution to the neighboring grid points
        Note: lots of redundant FLOPS here...should be optimized.*/
            for (l1=(1-M_sp);l1<=M_sp;l1++)/*loop over x dimension*/
            {   /*Note: We need to make sure the convolution wraps around 
                at the boundaries
                Calculate the boundary tests:*/
                lx=(m1+l1+M_rxd2)>=0;/*lx is true when reference is above 
                                     the lower x boundary*/
                rx=(m1+l1)<M_rxd2;   /*rx is true when reference is below 
                                     the upper x boundary*/
                xind = m1+l1+(rx-lx)*M_rx + M_rxd2;/*number in [0, M_rx -1]*/
                E_23x=E_2x[M_sp+l1-1]*E_3x[M_sp+l1-1];
                outArrayR[i]+= V0r*E_23x*ftaur[xind];
                outArrayI[i]+= V0r*E_23x*ftaui[xind];
            }
    }
    free(E_2x);
    return;
}