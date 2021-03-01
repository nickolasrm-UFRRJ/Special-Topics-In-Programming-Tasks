#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#define ERROR 1E-16
#define EPSILON 1E-17
struct stMatrix{
    double *A;
    double *LU;
    int mn;
};
typedef struct stMatrix tpMatrix;
void LoadMatrix(char *matrixFile, tpMatrix *matrix);
void PrintMatrix(const tpMatrix *matrix);
void LUMethod(tpMatrix *matrix);
void PrintX(double *X, const int size);


int main(int ac, char**av) {
  tpMatrix matrix;
  
  int      flagSave = atoi(av[3]);
  fprintf(stdout, "\nDecomposição LU\n");
  
  matrix.mn = atoi(av[2]);
  
  matrix.A = (double*) malloc (matrix.mn * matrix.mn *  sizeof(double));
  matrix.LU = (double*) malloc (matrix.mn * matrix.mn *  sizeof(double));
  memset(matrix.LU, 0x00, matrix.mn * matrix.mn *  sizeof(double));
  /*
  matrix.A[0]  = 2.0; matrix.A[1]  =  3.0; matrix.A[2]  =  1.0;  matrix.A[3]  =  5.0;
  matrix.A[4]  = 6.0; matrix.A[5]  = 13.0; matrix.A[6]  =  5.0;  matrix.A[7]  = 19.0;
  matrix.A[8]  = 2.0; matrix.A[9]  = 19.0; matrix.A[10] = 10.0;  matrix.A[11] = 23.0;
  matrix.A[12] = 4.0; matrix.A[13] = 10.0; matrix.A[14] = 11.0;  matrix.A[15] = 31.0;
  */
  
  LoadMatrix(av[1], &matrix);
 //PrintLU(matrix.A, matrix.mn);
  LUMethod(&matrix);
  //LUckechk(matrix);
  if (flagSave == 1)
    PrintX(matrix.LU, matrix.mn);


  free(matrix.LU);
  
  return EXIT_SUCCESS;
}

/*
 * Carrega a matrix e o vetor B do arquivo
 */
void LoadMatrix(char *matrixFile,tpMatrix *matrix){
  FILE *ptr = fopen(matrixFile, "rb+");
  assert(ptr != NULL);
  fread (matrix->A,sizeof(double), matrix->mn * matrix->mn, ptr);
  fclose(ptr);

  
}
/*
void PrintMatrix(const tpMatrix *matrix){

  fprintf(stdout, "Matrix (%d, %d)\n", matrix->mn, matrix->mn);
  for (int j = 0; j < matrix->mn; j++){
    for (int i = 0; i < matrix->mn; i++){
      int k = j * matrix->mn + i;
      fprintf(stdout, "%.7f ", matrix->LU[k]);
    }
    fprintf(stdout, "\n");
  }
}
*/

void calculate_st_end(int *stvar, int *endvar, int st, int end)
{
  int id = omp_get_thread_num();
  int nthreads = omp_get_num_threads();
  int interval = end - st;
  stvar[0] = id * interval /nthreads + st;
  endvar[0] = (id+1) * interval / nthreads + st;
}

void LUMethod(tpMatrix *matrix){
  int r = matrix->mn,
      c = matrix->mn;
      
  double *LU = matrix->LU,
         *A  = matrix->A,
          a = 0.0,
          u = 0.0,
          l = 0.0;
         
  l = matrix->A[0];
  //caso especial 1a linha e 1a coluna
  LU[0] = A[0];
  for (int i = 1; i < c; i++){
    LU[i] = A[i];
    LU[i * c] = A[i * c] / l;
  }//end-for (int i = 0; i < c; i++){
  //-----------------------------------------------------------------------
  
  for (int j = 1; j < r; j++){
    //calcular matrix U
    #pragma omp parallel firstprivate(LU, A) private(u, a)
    {
      int st, end;
      calculate_st_end(&st, &end, j, c);
      double *Ajc = A + j * c;
      double *LUjc = LU + j * c;
      double *LU_i = LU + st;
      double *LUkc_i;

      for (int i = st; i < end; i++, LU_i++){
        a = Ajc[i]; //i = j
        u = 0.0;
        LUkc_i = LU_i;
        for (int k = 0; k < j; k++, LUkc_i += c){
            u += LUjc[k] * LUkc_i[0];
        }
        u = a - u;
        LUjc[i] = u;
      }
    }
    
     //PrintMatrix(matrix);
    //calcular matrix L
     u = LU[j * c + j];
    #pragma omp parallel firstprivate(LU, A) private(l, a)
    {
      int st, end;
      calculate_st_end(&st, &end, j+1, r);

      double *LUkc = LU + st*c;
      double *Akc = A + st*c;
      double *LUic = LU;
      for (int k = st; k < end; k++, LUkc += c, Akc += c){
          l = 0.0;
          for (int i = 0; i < j; i++, LUic += c){
              l += LUkc[i] * LUic[j];
          }
          LUic = LU;
          a = Akc[j];
          l = (a - l) / u;
          LUkc[j] = l;
      }
    }
  //  PrintMatrix(matrix);
      
  }//end-for (int j = 1; j < r; j++{
  
  
  
}

void PrintX(double *X, const int size){
  int mn = size * size;
  FILE *ptr = fopen("solucao.bin", "w+");
  assert(ptr != NULL);
  fwrite (&mn , sizeof(int), 1, ptr);
  fwrite (X , sizeof(double), mn, ptr);
  fclose(ptr);
}
