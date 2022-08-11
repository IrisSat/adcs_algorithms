//gcc -o main.exe main.c matrix_math.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_math.h"


int main(int argc, char ** argv)
{
    /*** Matrix math ***/
    /*
    int m = 2;
    int n = 3;
    int p = 2;
    float A[2][3] = {
        {1,2,3},
        {4,5,6}
    };
    float B[3][2] = {
        {7,8},
        {9,10},
        {11,12}
    };
    printf("*** Matrix math ***\n\n");
    printf("Matrix A:\n");
    print_matrix(m,n,A);
    printf("Matrix B:\n");
    print_matrix(n,p,B);

    // Matrix multiplication example
    float C[2][2];
    matmul(m,n,p,A,B,C);
    printf("\nAxB:\n");
    print_matrix(n,p,C);
    printf("\n\n");

    // Transpose example
    float A_T[3][2];
    transpose(m,n,A,A_T);
    printf("A transpose:\n");
    print_matrix(n,m,A_T);
    printf("\n\n");

    //Vector Math
    float a[3] = {2,3,4};
    float b[3] = {5,6,7};
    printf("*** Vector math ***\n");
    printf("a = ");
    print_vector(3,a);
    printf("b = ");
    print_vector(3,b);

    // Dot product example
    float c_dot;
    dot_prod(a,b,&c_dot);
    printf("\na . b = %.2f\n",c_dot);

    // Cross product example
    float c_cross[3];
    cross_prod(a,b,c_cross);
    printf("a x b = ");
    print_vector(3,c_cross);
    printf("\n");
    */

    // Linear system solver example
    float A2[5][5] = {
        {4, 9, 3, 5, 8},
        {2, 1, 8, 0, 2},
        {3, 1, 4, 1, 5},
        {9, 2, 6, 5, 3},
        {5, 8, 9, 7, 9}
    };
    float b2[5] = {2, 3, 4, 5, 6}; // Technically this should be a column vector for the dims to allign but solve assumes it is already
    float x2[5]; //soln
    for(long solves = 0; solves<100000000; solves++){
        solve(5, A2, b2, x2);
        //solveInv(5, A2, b2, x2);
    }
    return 0;
}