// ADCS algorithms Version 4
#include <stdio.h>
#include <math.h>
void Eulers(float *I1, float *I2, float *I3, float *omega, float *T, float *omega_dot );


void EKF();

void Calc_Jacobian(float *I1, float *I2, float *I3, float *height, float *omega, float *S_p, float Jacobian_A[][5]);
void CalcQd(float *I1, float *I2, float *I3, float *w, float Qd[][5]);
void PredictP(float Jacobian_A[][5], float Pk[][5], float Qd[][5]);
void CalcK(float Pk[][5], float R[][5], float K[][5]); 

//***************************Matrix Math Functions***********************// 
void matmul(int m, int n, int p, float A[m][n], float B[n][p], float C[m][p]);
void transpose(int m, int n, float A[][n], float A_T[][m]);
void addmatrix(int m, int n, float A[][n], float B[][n], float C[][n]);
void vectorMultiply(int m, int n, float A[m][n], float B[m], float C[m]);
void vector_sub(int m, float *a, float *b, float *c);
void cross_prod(float *a, float *b, float *c);
void vector_mult_scalar(float *v, float mult);
void vector_add(int m, float *a, float *b, float *c);
void invert(int m, float A[m][m], float AInv[m][m]);
void solveInv(int m, float A[m][m], float B[m][m], float X[m][m]);
//************************************************************************//

// Variables required for Euler's Equations 
float I1=3.6287e-2; float I2=3.4092e-2; float I3=0.9026e-2; 
float omega_dot[3];
float T[3] ={0,0,0};
float omega[3]={0.01745, 0.001745, 0.001745}; 



int main(){
EKF();	
}


void EKF(){
// Variables required for EKF 
float w=1e-5; //Assumed process noise 
float height=0.01; 

float Pk[5][5]={
 {0.0076154,0,0,0,0}, 
 {0,0.0076154,0,0,0}, 
 {0,0,0.0076154,0,0}, 
 {0,0,0,304.617,0}, 
 {0,0,0,0,304.617}
};

float H[5][5]={
{1,0,0,0,0},
{0,1,0,0,0},
{0,0,1,0,0},
{0,0,0,1,0},
{0,0,0,0,1}
};

float R[5][5]={
{7.61543E-5,0,0,0,0},
{0,7.61543E-5,0,0,0},
{0,0,7.61543E-5,0,0},
{0,0,0,7.61543E-5,0},
{0,0,0,0,7.61543E-5}
};

float S_p[2]={0,0};
float Jacobian_A[5][5]={
{0,0,0,0,0},
{0,0,0,0,0},
{0,0,0,0,0},
{0,0,0,0,0},
{0,0,0,0,0}
};

float Qd[5][5];	
float K [5][5];


Calc_Jacobian(&I1, &I2, &I3, &height, omega, S_p, Jacobian_A);
CalcQd(&I1, &I2, &I3, &w, Qd);
PredictP(Jacobian_A,  Pk, Qd);
printf("Pk is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[0][0], Pk[0][1], Pk[0][2], Pk[0][3],Pk[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[1][0], Pk[1][1], Pk[1][2], Pk[1][3],Pk[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[2][0], Pk[2][1], Pk[2][2], Pk[2][3],Pk[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[3][0], Pk[3][1], Pk[3][2], Pk[3][3],Pk[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[4][0], Pk[4][1], Pk[4][2], Pk[4][3],Pk[4][4]);
CalcK(Pk, R, K);

	
}


void Eulers(float *I1, float *I2, float *I3, float *omega, float *T, float *omega_dot ){
	omega_dot[0]=(T[0]-(*I3-*I2)*omega[1]*omega[2])/(*I1);
	omega_dot[1]=(T[1]-(*I1-*I3)*omega[0]*omega[2])/(*I2);
    omega_dot[2]=(T[2]-(*I2-*I1)*omega[0]*omega[1])/(*I3); 
}



void Calc_Jacobian(float *I1, float *I2, float *I3, float *height, float *omega, float *S_p, float Jacobian_A[][5]){
float xp=S_p[0]; float yp=S_p[1]; 
float w1=omega[0]; float w2=omega[1]; float w3=omega[2]; 

Jacobian_A[0][1]=-((*I3-*I2)/(*I1))*w3;
Jacobian_A[0][2]=-((*I3-*I2)/(*I1))*w2;

Jacobian_A[1][0]=-((*I1-*I3)/(*I2))*w3; 
Jacobian_A[1][2]=-((*I1-*I3)/(*I2))*w1; 

Jacobian_A[2][0]=-((*I2-*I1)/(*I3))*w2;
Jacobian_A[2][1]=-((*I2-*I1)/(*I3))*w1;

Jacobian_A[3][1]=*height; 
Jacobian_A[3][4]=-w3;

Jacobian_A[4][0]=-*height; 
Jacobian_A[4][3]=w3; 
//Jacobian_A[5][5]={
//{0,   -((*I3-*I2)/(*I1))*w3,  -((*I3-*I2)/(*I1))*w2, 0, 0 },
//{-((*I1-*I3)/(*I2))*w3,  0,   -((*I1-*I3)/(*I2))*w1, 0, 0 },
//{-((*I2-*I1)/(*I3))*w2, -((*I2-*I1)/(*I3))*w1, 0,    0, 0 } ,
//{-yp*xp/(*height),  (*height)+(xp*xp)/(*height),  -yp,  (-yp*w1/(*height))+(2*xp*w2/(*height)), -w3-(xp*w1/(*height))},
//{-(*height)-(yp*yp/(*height)),     xp*yp/(*height),  xp,  w3+(yp*w2/(*height)), (2*yp*w1)/(*height)+(xp*w2/(*height))}
//};
}


void CalcQd(float *I1, float *I2, float *I3, float *w, float Qd[][5]){
	float Bd[5][3]={
	{1/(*I1), 0, 0},
	{0, 1/(*I2), 0}, 
	{0, 0, 1/(*I3)},
	{0, 0, 0, },
    {0, 0, 0, }
	};
	
	float W[3][3]={
	{*w*(*w),0,0},
	{0,*w*(*w),0},
	{0,0,*w*(*w)}
	};
	
	float C[5][3];
	float Bd_Transpose[3][5];
	matmul(5,3,3,Bd, W, C);
	transpose(5,3,Bd, Bd_Transpose);
	matmul(5,3,5,C,Bd_Transpose, Qd);	
}


void PredictP(float Jacobian_A[][5], float Pk[][5], float Qd[][5]){
	// This function calculates & updates the error covariance (Pk)
	float C[5][5];
	float D[5][5];
    float E[5][5];
    
	matmul(5,5,5,Jacobian_A, Pk, C);
	transpose(5,5, Jacobian_A,D);
	matmul(5,5,5,C,D,E);
	addmatrix(5,5,E,Qd,Pk);
	printf("Pk is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[0][0], Pk[0][1], Pk[0][2], Pk[0][3],Pk[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[1][0], Pk[1][1], Pk[1][2], Pk[1][3],Pk[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[2][0], Pk[2][1], Pk[2][2], Pk[2][3],Pk[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[3][0], Pk[3][1], Pk[3][2], Pk[3][3],Pk[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[4][0], Pk[4][1], Pk[4][2], Pk[4][3],Pk[4][4]);
}

void CalcK(float Pk[][5], float R[][5], float K[][5]){
	
//The equation to calculate K simplifies to Pk/(Pk+R) if H is just the identity matrix. 
float Pk_add_R[5][5]; 

addmatrix(5, 5, Pk, R, Pk_add_R );
solveInv(5, Pk, Pk_add_R, K);

printf("Pk is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[0][0], Pk[0][1], Pk[0][2], Pk[0][3],Pk[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[1][0], Pk[1][1], Pk[1][2], Pk[1][3],Pk[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[2][0], Pk[2][1], Pk[2][2], Pk[2][3],Pk[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[3][0], Pk[3][1], Pk[3][2], Pk[3][3],Pk[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[4][0], Pk[4][1], Pk[4][2], Pk[4][3],Pk[4][4]);

printf("K is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", K[0][0], K[0][1], K[0][2], K[0][3],K[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[1][0], K[1][1], K[1][2], K[1][3],K[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[2][0], K[2][1], K[2][2], K[2][3],K[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[3][0], K[3][1], K[3][2], K[3][3],K[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", K[4][0], K[4][1], K[4][2], K[4][3],K[4][4]);

	
}




//***************************Matrix Math Functions***********************// 

void vector_mult_scalar(float *v, float mult)
{
    for(int i=0; i < 3; i++)
        v[i] *= mult;
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

void transpose(int m, int n, float A[][n], float A_T[][m])
{
    int i,j;
    for(i=0; i < m; i++){
        for(j=0; j < n; j++){
            A_T[j][i] = A[i][j];
        }
    }
}

void addmatrix(int m, int n, float A[][n], float B[][n], float C[][n])
{
    int i,j;
    for(i=0; i < m; i++){
        for(j=0; j < n; j++){
            C[j][i] = A[j][i]+B[j][i];
        }
    }
}

void vectorMultiply(int m, int n, float A[m][n], float B[m], float C[m])
{
    int i,j;
    for(i=0; i < m; i++){
        for(j=0; j < n; j++){
            C[i] += A[i][j]*B[j];
            
        }
    }
}

void vector_sub(int m, float *a, float *b, float *c)
{
    for(int i=0; i < m; i++)
        c[i] = a[i]-b[i];
}

void cross_prod(float *a, float *b, float *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
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

void solveInv(int m, float A[m][m], float B[m][m], float X[m][m]){
    // This solves the equation XA = B which is equivalent to X = B / A in matlab
    float Ainv[m][m];
    invert(m, A, Ainv);
    // Ax = b === x = A^-1 * b
    
    
    matmul(m, m, m, B, Ainv, X);
    //print_matrix(m, 1, xmat);
    
}