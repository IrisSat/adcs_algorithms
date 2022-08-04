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

void solve(int m, float A[m][m], float b[m], float x[m]){
    float L[m][m];
    float U[m][m];
    float Y[m];
    int n = m; 
    int j = 0;
    int i = 0;
    int k = 0;
    
    // set all values to zero
    for(int i=0; i < m; i++){
        for(int j=0; j < m; j++){
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }

    for(i=0; i<n; i++){
        //upper
        for(k = i; k<n; k++){
            float sum = 0;
            for(j=0; j < i; j++){
                sum += L[i][j] * U[j][k];
            }
            U[i][k] = A[i][k] - sum;
        }
        //lower
        for(k = i; k<n; k++){
            if(i==k){
                L[i][i] = 1;
            }
            else{
                float sum = 0;
                for(j = 0; j< i; j++){
                    sum += L[k][j] * U[j][i];
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }

    for(i=0; i<n; i++)
    {
        Y[i]=b[i];
        for(j=0; j<i; j++)
        {
            Y[i]-=L[i][j]*Y[j];
        }
    }

    for(i=n-1; i>=0; i--)
    {
        x[i]= Y[i];
        for(j=i+1; j<n; j++)
        {
            x[i]-=U[i][j]*x[j];
        }
        x[i]/=U[i][i];
    }
    //print_vector(m, x);
    /*verify soln
    for(i = 0; i < n; i++){
        printf("Checking equation %i\n", i);
        float temp = 0;
        for(j = 0; j < n; j++){
            temp += A[i][j] * x[j];
        }
        if(temp != b[i]){
            printf("Welp, this one failed since %f != %f\n", temp, b[i]);
        }
    }*/
}
void invert(int m, float A[m][m], float AInv[m][m]){
    int i,j,k,matsize;
    matsize = m;
    float temp;
    for(i=0;i<matsize;i++){                               //automatically initialize the unit matrix, e.g.
        for(j=0;j<matsize;j++){                            //  -       -
            if(i==j)                                        // | 1  0  0 |
                AInv[i][j]=1;                                  // | 0  1  0 |
            else                                            // | 0  0  1 |
                AInv[i][j]=0; 
        }
    }
    for(k=0;k<matsize;k++)                                  //by some row operations,and the same row operations of
    {                                                       //Unit mat. I gives the inverse of matrix A
        temp=A[k][k];                   //'temp'  
        // stores the A[k][k] value so that A[k][k]  will not change
        for(j=0;j<matsize;j++)      //during the operation //A[i] //[j]/=A[k][k]  when i=j=k
        {
            A[k][j]/=temp;                                  //it performs // the following row operations to make A to unit matrix
            AInv[k][j]/=temp;                                  //R0=R0/A[0][0],similarly for I also
                                                            //R0=R0/A[0][0]
        }                                                   //R1=R1-R0*A[1][0] similarly for I
        for(i=0;i<matsize;i++)                              //R2=R2-R0*A[2][0]      ,,
        {
            temp=A[i][k];                       //R1=R1/A[1][1]
            for(j=0;j<matsize;j++)             //R0=R0-R1*A[0][1]
            {                                   //R2=R2-R1*A[2][1]
                if(i==k)
                    break;                      //R2=R2/A[2][2]
                A[i][j]-=A[k][j]*temp;          //R0=R0-R2*A[0][2]
                AInv[i][j]-=AInv[k][j]*temp;          //R1=R1-R2*A[1][2]
            }
        }
    }
}

void solveInv(int m, float A[m][m], float b[m], float x[m]){
    float Ainv[m][m];
    invert(m, A, Ainv);
    // Ax = b === x = A^-1 * b
    float xmat[m][1];
    float bmat[m][1];
    for(int i=0; i<m; i++){
        bmat[i][0] = b[i];
    }
    
    matmul(m, m, 1, Ainv, bmat, xmat);
    //print_matrix(m, 1, xmat);
    for(int i = 0; i < m; i++){
        x[i] = xmat[i][0];
    }
}

