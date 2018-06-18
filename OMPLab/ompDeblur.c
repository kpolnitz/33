//OpenMP version.  Edit and submit only this file.
/* Enter your details below
 * Name : Kaela Polnitz
 * UCLA ID: 504751829
 * Email id: kaelanpolnitz@gmail.com
 Input - Old Files 
 */
//ROUGH SPEED: 11.148525x speed-up
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

int OMP_xMax;
#define xMax OMP_xMax
int OMP_yMax;
#define yMax OMP_yMax
int OMP_zMax;
#define zMax OMP_zMax

int OMP_Index(int x, int y, int z)
{
	return ((z * yMax + y) * xMax + x);
}
#define Index(x, y, z) OMP_Index(x, y, z)

double OMP_SQR(double x)
{
	return pow(x, 2.0);
}
#define SQR(x) OMP_SQR(x)

double* OMP_conv;
double* OMP_g;

void OMP_Initialize(int xM, int yM, int zM)
{
	xMax = xM;
	yMax = yM;
	zMax = zM;
	assert(OMP_conv = (double*)malloc(sizeof(double) * xMax * yMax * zMax));
	assert(OMP_g = (double*)malloc(sizeof(double) * xMax * yMax * zMax));
}
void OMP_Finish()
{
	free(OMP_conv);
	free(OMP_g);
}
void OMP_GaussianBlur(double *u, double Ksigma, int stepCount)
{
    int l = 1;
    int xLim = xMax - 2;
    int yLim = yMax - 2;
    int zLim = zMax - 2;
	double lambda = (Ksigma * Ksigma) / (double)(stepCount << 1);
	double nu = (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
    int x, y, z, step, iz, iy;
	double boundryScale = 1.0 / (1.0 - nu);
	double postScale = pow(nu / lambda, (double)((stepCount << 1) + (stepCount << 0)));
    int gXY = xMax*yMax;
    int c = (zMax-1)*gXY;
    int g = (yMax-1)*xMax;
    int s = boundryScale*postScale;
    int t = zLim*gXY;
    int zz = yLim*xMax;
	//for(step = 0; step < stepCount; step+=4)
	//{

        iy = 0;

        iz = gXY;
        for(z = 1; z < zMax; z++)
        {
            
            for(x = 0; x < xMax; x++)
            {
                iy = 0;
                for(y = 0; y < yMax; y++)
                {
                    int o = iy + iz + x;
					u[o] = u[o - gXY] * nu;
                    iy += xMax;
				}
                
                
			}
            
            iy = 0;

            
            iz += gXY;
		}

        
	//} //end of step loop
    



}
void OMP_Deblur(double* u, const double* f, int maxIterations, double dt, double gamma, double sigma, double Ksigma)
{
    //let's use a flag
    int didfirstloop = 3;
	double epsilon = 1.0e-7;
	double sigma2 = SQR(sigma);
    int x, y, z, iteration, iz, iy;
	int converged = 0;
	int lastConverged = 0;
	int fullyConverged = (xMax - 1) * (yMax - 1) * (zMax - 1);
	double* conv = OMP_conv;
	double* g = OMP_g;
    int xymult = xMax*yMax;
	for(iteration = 0; iteration < maxIterations && converged != fullyConverged; iteration++)
    {   iz = xymult;
        
        for(z = 1; z < zMax - 1; z++)
        {
            
            for(x = 1; x < xMax - 1; x++)
            {
                iy = xMax;
                //y = 1;
                for(y = 1; y < yMax - 1; y++)
                {
                    
                    int izyx = iz + iy + x;
                    double access = u[izyx];
					g[izyx] = 1.0 / sqrt(epsilon +
						SQR(access - u[izyx+1]) +
						SQR(access - u[izyx-1]) +
						SQR(access - u[izyx+xMax]) +
						SQR(access - u[izyx-xMax]) +
						SQR(access - u[izyx+xymult]) +
						SQR(access - u[izyx-xymult]));
                    iy += xMax;
				}
			}
            iz += xymult;
            
		}
		memcpy(conv, u, sizeof(double) * xMax * yMax * zMax);
		OMP_GaussianBlur(conv, Ksigma, 3);

		//OMP_GaussianBlur(conv, Ksigma, 3);
		converged = 0;
        int izzz = xymult;
        
        for(z = 1; z < zMax - 1; z++)
        {
            
            int ab = 0;
            //for(x = 0; x < xMax; x++)
            //{
            iy = 0;
            //for(y = 0; y < yMax; y++)
            //{
            
            int izzyyx = iz + iy + ab;
            double r = conv[izzyyx] * f[izzyyx] / sigma2;
            r = (r * (2.38944 + r * (0.950037 + r))) / (4.65314 + r * (2.57541 + r * (1.48937 + r)));
            conv[izzyyx] -= f[izzyyx] * r;
            iy += xMax;
            //}
            //}
            ab++;
            iz += xymult;
            
            for(x = 1; x < xMax - 1; x++)
            {
                //iyyy can go here, and then iterate by += !!!
                int iyyy = xMax;
                for(y = 1; y < yMax - 1; y++)
                {
                    
                    int zystuff = izzz+iyyy;
                    int xyzz = zystuff+x;
                    int i1 = xyzz-1;
                    int i2 = i1+2;
                    int i3 = xyzz-xMax;
                    int i4 = xyzz+xMax;
                    int i5 = xyzz-xymult;
                    int i6 = xyzz+xymult;
					double oldVal = u[xyzz];
					double newVal = (u[xyzz] + dt * (
						u[i1] * g[i1] +
						u[i2] * g[i2] +
						u[i3] * g[i3] +
						u[i4] * g[i4] +
						u[i5] * g[i5] +
						u[i6] * g[i6] - gamma * conv[xyzz])) /
						(1.0 + dt * (g[xyzz+1] + g[xyzz-1] + g[i4] + g[i3] + g[i6] + g[i5]));
					if(fabs(oldVal - newVal) < epsilon)
					{
						converged++;
					}
					u[xyzz] = newVal;
                    iyyy+=xMax;
				}
			}
            izzz+=xymult;
		}
		if(converged > lastConverged)
		{
			printf("%d pixels have converged on iteration %d\n", converged, iteration);
			lastConverged = converged;
		}
	}
}

