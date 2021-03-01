#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#define STABILITY 1.0f/sqrt(3.0f)


void mdf_heat(double ***  __restrict__ u0,
              double ***  __restrict__ u1,
              const unsigned int npX,
              const unsigned int npY,
              const unsigned int npZ,
              const unsigned int ix,
              const unsigned int jy,
              const unsigned int kz,
              const double deltaH,
              const double deltaT,
              const unsigned int tsteps,
              const double boundaries);


void save2bin(double ***u,
             const unsigned int npX,
             const unsigned int npY,
             const unsigned int npZ);

int main (int ac, char **av){
  double ***u0;
  double ***u1;
  double deltaT = 0.01;
  double deltaH = 0.25f;
  
  unsigned int tsteps = (unsigned int) atoi(av[1]); //(time / deltaT);
  unsigned int npX = (unsigned int) atoi(av[2]); //(sizeX / deltaH); //Number of points in X axis
  unsigned int npY = (unsigned int) atoi(av[3]); //(sizeY / deltaH);
  unsigned int npZ = (unsigned int) atoi(av[4]); //(sizeZ / deltaH);


  double temp   = atof(av[5]);
  int flag2save = atoi(av[6]);

  fprintf(stdout, "\nSimulação - Domínio(x = %u, y = %u, z = %u, t = %u)\n", npX, npY, npZ, tsteps);

  u0 = (double***) malloc (npZ * sizeof(double**));
  u1 = (double***) malloc (npZ * sizeof(double**));

  for (unsigned int i = 0; i < npZ; i++){
    u0[i] = (double**) malloc (npY * sizeof(double*));
    u1[i] = (double**) malloc (npY * sizeof(double*));
  }


  for (unsigned int i = 0; i < npZ; i++){
    for (unsigned int j = 0; j < npY; j++){
      double *aux0 = (double *) malloc (npX * sizeof(double));
      double *aux1 = (double *) malloc (npX * sizeof(double));
      //initial condition - zero in all points
      memset(aux0, 0x00, npX * sizeof(double));
      memset(aux1, 0x00, npX * sizeof(double));
      u0[i][j] = aux0;
      u1[i][j] = aux1;
    }
  }


  mdf_heat(u0, u1, npX, npY, npZ, npX/2, npY/2, npZ/2, deltaH, deltaT, tsteps, temp);

  if (flag2save == 1)
    save2bin(u0, npX, npY, npZ);


  //Free memory
  for (unsigned int i = 0; i < npZ; i++){
    for (unsigned int j = 0; j < npY; j++){
      free(u0[i][j]);
      free(u1[i][j]);
    }
  }

  for (unsigned int i = 0; i < npZ; i++){
    free(u0[i]);
    free(u1[i]);
  }



  free(u0);
  free(u1);


  return EXIT_SUCCESS;
}



void mdf_heat(double ***  __restrict__ u0,
              double ***  __restrict__ u1,
              const unsigned int npX,
              const unsigned int npY,
              const unsigned int npZ,
              const unsigned int ix,
              const unsigned int jy,
              const unsigned int kz,
              const double deltaH,
              const double deltaT,
              const unsigned int tsteps,
              const double boundaries){

    const double alpha = deltaT / (deltaH * deltaH);
    assert(alpha < STABILITY);

    unsigned int step = tsteps / 20;
    if (step == 0) step = 1;

    for ( unsigned int steps = 0; steps < tsteps; steps++){
              #pragma omp parallel for
              for (unsigned int i = 0; i < npZ; i++){
                for (unsigned int j = 0; j < npY; j++){
                  for (unsigned int k = 0; k < npX; k++){
                      register double left   = 0.0;
                      register double right  = 0.0;
                      register double up     = 0.0;
                      register double down   = 0.0;
                      register double top    = 0.0;
                      register double bottom = 0.0;
                      register double center = u0[i][j][k];
                      
                      if ((i == ix) && (j == jy) && (k == kz)){
                        center = boundaries;
                        u0[i][j][k] = boundaries;
                      }
                        

                      if ((k > 0) && (k < (npX - 1))){
                        left  = u0[i][j][k-1];
                        right = u0[i][j][k+1];
                      }else if (k == 0){
                        right = u0[i][j][k+1];
                        left = right;
                      }else{
                        left = u0[i][j][k-1];
                        right = left;
                        
                      } 

                      if ((j > 0) && (j < (npY - 1))){
                        up  = u0[i][j-1][k];
                        down = u0[i][j+1][k];
                      }else if (j == 0){
                        down = u0[i][j+1][k];
                        up = down;
                      } 
                      else {
                        up = u0[i][j-1][k];
                        down = up;
                      }
                        

                      if ((i > 0) && (i < (npZ - 1))){
                        top  = u0[i-1][j][k];
                        bottom = u0[i+1][j][k];
                      }else if (i == 0){
                        bottom = u0[i+1][j][k];
                        top = bottom;
                      } 
                      else{
                        top = u0[i-1][j][k];
                        bottom = top;
                      } 


                      u1[i][j][k] =  alpha * ( top + bottom + up + down + left + right  - (6.0f * center )) + center;
                  }
                }
              }

              double ***ptr = u0;
              u0 = u1;
              u1 = ptr;

              if ((steps % step) == 0){
                fprintf(stdout, ".");
                fflush(stdout);
              }

    }
}

void save2bin(double ***u,
             const unsigned int npX,
             const unsigned int npY,
             const unsigned int npZ){
  char fileName[128];
  sprintf(fileName, "%s-%d-%d-%d-log.bin", __FILE__,  npX, npY, npZ);
  fprintf(stdout, "\nSaving file [%s] ", fileName); fflush(stdout);
  FILE *ptr = fopen(fileName, "wb+");
  assert(ptr != NULL);

   for (unsigned int i = 0; i < npZ; i++){
     for (unsigned int j = 0; j < npY; j++){
         fwrite ((const void*)u[i][j] , sizeof(double), npX, ptr);
     }
   }

   fprintf(stdout, "\t[OK]\n");
   fclose(ptr);

}

