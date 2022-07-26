//gcc -o test.exe test_matrix_math.c matrix_math.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "matrix_math.h"

#define MAX_LINE_LENGTH 80

float eps = 1e-3;

bool float_compare(float a, float b);
int test_matmul(bool verbose);
int test_cross_prod(bool verbose);

int main(int argc, char ** argv)
{
    test_matmul(true);
    test_cross_prod(true);
    return 0;
}


int test_matmul(bool verbose)
{
    printf("Testing matmul:\n\n");
    FILE *fp;

    char buf[MAX_LINE_LENGTH] = {0};
    // int num_tests;

    // Open matmul_test.txt
    char *path = "./test_files/matmul_test.txt";
    fp = fopen(path, "r");
    if(!fp){
        perror(path);
        return 1;
    }
    // Clear first line
    fgets(buf,MAX_LINE_LENGTH,fp);
    // Get number of tests
    int num_tests;
    if(!fscanf(fp,"%d",&num_tests)){
        printf("Error reading num_tests\n");
        return 1;
    }
    printf("num_tests: %d\n",num_tests);


    int i,j,n_test;
    int m,n,p;
    for(n_test=0; n_test < num_tests; n_test++){
        printf("Test %d:\n",n_test+1);
        // printf("buf: %s\n",buf);
        fgets(buf,MAX_LINE_LENGTH,fp); // newline
        fgets(buf,MAX_LINE_LENGTH,fp); // Test number
        fgets(buf,MAX_LINE_LENGTH,fp); // m n p
        if(!fscanf(fp,"%d",&m)){
            printf("Error reading m\n");
            return 1;
        }
        if(!fscanf(fp,"%d",&n)){
            printf("Error reading n\n");
            return 1;
        }
        if(!fscanf(fp,"%d",&p)){
            printf("Error reading p\n");
            return 1;
        }
        printf("m: %d, n: %d, p: %d\n",m,n,p);
        fgets(buf,MAX_LINE_LENGTH,fp); // newline after m n p
        fgets(buf,MAX_LINE_LENGTH,fp); // A
        // Initialize matrices
        float A[m][n];
        for(i=0; i < m; i++){
            for(j=0; j < n; j++){
                // fscanf(fp,"%f",&A[i][j]);
                if(!fscanf(fp,"%f",&A[i][j])){
                    printf("Error reading matrix A\n");
                    return 1;
                }
            }
        }
        fgets(buf,MAX_LINE_LENGTH,fp); // newline after final row of A
        fgets(buf,MAX_LINE_LENGTH,fp); // B
        float B[n][p];
        for(i=0; i < n; i++){
            for(j=0; j < p; j++){
                // fscanf(fp,"%f",&B[i][j]);
                if(!fscanf(fp,"%f",&B[i][j])){
                    printf("Error reading matrix B\n");
                    return 1;
                }
            }
        }
        fgets(buf,MAX_LINE_LENGTH,fp); // newline after final row of B
        fgets(buf,MAX_LINE_LENGTH,fp); // C
        float C[m][p];
        for(i=0; i < m; i++){
            for(j=0; j < p; j++){
                // fscanf(fp,"%f",&C[i][j]);
                if(!fscanf(fp,"%f",&C[i][j])){
                    printf("Error reading matrix C\n");
                    return 1;
                }
            }
        }
        // Perform calculation
        float C_result[m][p];
        matmul(m,n,p,A,B,C_result);
        // Compare
        bool pass = true;
        for(i=0; i < m; i++){
            for(j=0; j < p; j++){
                if(!float_compare(C_result[i][j],C[i][j])){
                    pass = false;
                }
            }
        }
        if(verbose){
            printf("Matrix A:\n");
            print_matrix(m,n,A);
            printf("Matrix B:\n");
            print_matrix(n,p,B);
            printf("Matrix C (actual):\n");
            print_matrix(m,p,C);
            printf("Matrix C (calculated):\n");
            print_matrix(m,p,C_result);
        }
        printf("***Test %d pass: %d\n\n",n_test+1,pass);
    }

    if(fclose(fp)){
        perror(path);
        return 1;
    }

    return 0;
}

int test_cross_prod(bool verbose)
{
    printf("Testing cross_prod:\n\n");
    FILE *fp;

    char buf[MAX_LINE_LENGTH] = {0};
    // int num_tests;

    // Open matmul_test.txt
    char *path = "./test_files/crossprod_test.txt";
    fp = fopen(path, "r");
    if(!fp){
        perror(path);
        return 1;
    }
    // Clear first line
    fgets(buf,MAX_LINE_LENGTH,fp);
    // Get number of tests
    int num_tests;
    if(!fscanf(fp,"%d",&num_tests)){
        printf("Error reading num_tests\n");
        return 1;
    }
    printf("num_tests: %d\n",num_tests);


    int i,n_test;
    int n_dim = 3;
    for(n_test=0; n_test < num_tests; n_test++){
        printf("Test %d:\n",n_test+1);
        // printf("buf: %s\n",buf);
        fgets(buf,MAX_LINE_LENGTH,fp); // newline
        fgets(buf,MAX_LINE_LENGTH,fp); // Test number
        fgets(buf,MAX_LINE_LENGTH,fp); // a
        // Initialize vectors
        float a[n_dim];
        for(i=0; i < n_dim; i++){
            if(!fscanf(fp,"%f",&a[i])){
                printf("Error reading vector a\n");
                return 1;
            }
        }
        fgets(buf,MAX_LINE_LENGTH,fp); // newline after a
        fgets(buf,MAX_LINE_LENGTH,fp); // b
        float b[n_dim];
        for(i=0; i < n_dim; i++){
            if(!fscanf(fp,"%f",&b[i])){
                printf("Error reading vector b\n");
                return 1;
            }
        }
        fgets(buf,MAX_LINE_LENGTH,fp); // newline after b
        fgets(buf,MAX_LINE_LENGTH,fp); // c
        float c[n_dim];
        for(i=0; i < n_dim; i++){
            if(!fscanf(fp,"%f",&c[i])){
                printf("Error reading vector c\n");
                return 1;
            }
        }
        // Perform calculation
        float c_result[n_dim];
        cross_prod(a,b,c_result);
        // Compare
        bool pass = true;
        for(i=0; i < n_dim; i++){
            if(!float_compare(c_result[i],c[i])){
                pass = false;
            }
        }
        if(verbose){
            printf("Vector a:\n");
            print_vector(n_dim,a);
            printf("Vector b:\n");
            print_vector(n_dim,b);
            printf("Vector c (actual):\n");
            print_vector(n_dim,c);
            printf("Vector c (calculated):\n");
            print_vector(n_dim,c_result);
        }
        printf("***Test %d pass: %d\n\n",n_test+1,pass);
    }

    if(fclose(fp)){
        perror(path);
        return 1;
    }

    return 0;
}


bool float_compare(float val, float compare)
{
     return compare - eps < val && val < compare + eps;
}