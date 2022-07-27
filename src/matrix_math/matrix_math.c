#include <stdio.h>

#include "matrix_math.h"

/*** Matrix Math ***/
void matmul_scalar(int m, int n, float A[][n], float mult)
{
    int i,j;
    for(i=0; i < m; i++){
        for(j=0; j < n; j++){
            A[i][j] *= mult;
        }
    }
}

void matdiv_scalar(int m, int n, float A[][n], float div)
{
    int i,j;
    for(i=0; i < m; i++){
        for(j=0; j < n; j++){
            A[i][j] /= div;
        }
    }
}

void transpose(int m, int n, float A[][n], float A_T[][m])
{
    int i,j;
    for(i=0; i < m; i++){
        for(j=0; j < n; j++){
            A_T[j][i] = A[i][j];
        }
    }
}

void matmul(int m, int n, int p, float A[m][n], float B[n][p], float C[m][p])
{
    int i,j,k;
    for(i=0; i < m; i++){
        for(j=0; j < p; j++){
            C[i][j] = 0;
            for(k=0; k < n; k++){
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}


/*** Vector Math  ***/
void vector_add(float *a, float *b, float *c)
{
    for(int i=0; i < 3; i++)
        c[i] = a[i]+b[i];
}
void vector_sub(float *a, float *b, float *c)
{
    for(int i=0; i < 3; i++)
        c[i] = a[i]-b[i];
}

void vector_mult_scalar(float *v, float mult)
{
    for(int i=0; i < 3; i++)
        v[i] *= mult;
}

void vector_div_scalar(float *v, float div)
{
    for(int i=0; i < 3; i++)
        v[i] /= div;
}

void dot_prod(float *a, float *b, float *c)
{
    *c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void cross_prod(float *a, float *b, float *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}


/*** Utilities ***/
void print_vector(int n, float *v)
{
    int i;
    for(i=0; i < n; i++)
        printf("%.2f\t",v[i]);
    printf("\n");
}
void print_matrix(int m, int n, float A[][n])
{
    int i,j;
    for(i=0; i < m; i++){
        for(j=0; j < n; j++){
            printf("%.2f\t",A[i][j]);
        }
        printf("\n");
    }
}