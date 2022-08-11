#include <stdio.h>
#include <math.h>

void Eulers(float *I1, float *I2, float *I3, float *omega, float *T, float *omega_dot );
void RK4();
void sun_tranformation(); 
void EKF();

void Calc_Jacobian(float *I1, float *I2, float *I3, float *height, float *omega, float *S_p);
void CalcQd(float *I1, float *I2, float *I3, float *w);
void PredictP();
void CalcK();
void controller(float *EKF_bodyrates, float *EKF_pointing);


//***************************Matrix Math Functions***********************// 
void matmul(int m, int n, int p, float A[m][n], float B[n][p], float C[m][p]);
void transpose(int m, int n, float A[][n], float A_T[][m]);
void addmatrix(int m, int n, float A[][n], float B[][n], float C[][n]);
void vectorMultiply(int m, int n, float A[m][n], float B[m], float C[m]);
void vector_sub(int m, float *a, float *b, float *c);
void cross_prod(float *a, float *b, float *c);
//vector_mult_scalar(float *v, float mult); 
void vector_add(int m, float *a, float *b, float *c);

void vector_mult_scalar(float *v, float mult)
{
    for(int i=0; i < 3; i++)
        v[i] *= mult;
}
//************************************************************************//


// Variables required for Euler's Equations 
float I1=3.6287e-2; float I2=3.4092e-2; float I3=0.9026e-2; 
float omega_dot[3];
float T[3] ={0,0,0};
float omega[3]={0.01745, 0.001745, 0.001745}; 
float h=0.1; //Step size (s)
float height=0.01; 

// Variables required for EKF 
float w=1e-5; //Assumed process noise 
//Make this into just actual numbers. 
//float Pk[5][5]={
// {0.0076154,0,0,0,0}, 
// {0,0.0076154,0,0,0}, 
// {0,0,0.0076154,0,0}, 
// {0,0,0,304.617,0}, 
// {0,0,0,0,304.617}
//};

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
int eclipse_flag=0; 
//when eclipse flag is 1, the spacecraft is in eclipse. 

int main(){

EKF();
return 0; 
}


void Calc_Jacobian(float *I1, float *I2, float *I3, float *height, float *omega, float *S_p){
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


void CalcQd(float *I1, float *I2, float *I3, float *w){
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

void PredictP(){
	// This function calculates & updates the error covariance (Pk)
	float C[5][5];
	float D[5][5];
    float E[5][5];
    
	matmul(5,5,5,Jacobian_A, Pk, C);
	transpose(5,5, Jacobian_A,D);
	matmul(5,5,5,C,D,E);
	addmatrix(5,5,E,Qd,Pk);
}

void CalcK(){
	
//The equation to calculate K simplifies to Pk/(Pk+R) if H is just the identity matrix. 
float Pk_add_R[5][5]; 
addmatrix(5, 5, Pk, R, Pk_add_R );

	
}


void Eulers(float *I1, float *I2, float *I3, float *omega, float *T, float *omega_dot ){
	omega_dot[0]=(T[0]-(*I3-*I2)*omega[1]*omega[2])/(*I1);
	omega_dot[1]=(T[1]-(*I1-*I3)*omega[0]*omega[2])/(*I2);
    omega_dot[2]=(T[2]-(*I2-*I1)*omega[0]*omega[1])/(*I3); 
//	printf("omegax_dot is \t omegay_dot is \t omegaz_dot is \n ");
//	printf("%E \t %E \t  %E", omega_dot[0], omega_dot[1], omega_dot[2]);
}


void EKF(){	
Calc_Jacobian(&I1, &I2, &I3, &height, omega, S_p);
CalcQd(&I1, &I2, &I3, &w);
PredictP();
CalcK(); 
	
float True_state[5]={omega[0],omega[1],omega[2],S_p[0],S_p[1]} ;
float Measured_state[5]={0.0221420309324296, 0.0177486658038343, -0.0179671575114262, 0.369903257448045, -0.380396852672751};


// Check for eclipse 
if (eclipse_flag ==1){
	H[3][3]=0;
	H[4][4]=0; 
}

//State updates 
float Xk[5]; 
float H_times_True_state[5];
vectorMultiply(5, 5, H, True_state, H_times_True_state);
float M_minus_T[5]; 
vector_sub(5, Measured_state, H_times_True_state, M_minus_T);

printf("%E \t %E \t %E \t %E \t %E \n  ", M_minus_T[0], M_minus_T[1], M_minus_T[2], M_minus_T[3],M_minus_T[4]);

//Multiply K*M_minus_T.
// Add True_state to product of previous step. 
}

void controller(float *EKF_bodyrates, float *EKF_pointing){
	//The controller takes in the estimated body rates and sun measurments determined by the EKF.
	//The controller outputs a command torque
	float desired_direction[3]={0,0,-1}; 
	float desired_rate[3]={0.01745329, 0, 0}; 
	float rate_gain=0.0004; float direction_gain=-0.000005; 
	float direction_command[3]; float rate_command[3]; 
	cross_prod(EKF_pointing, desired_direction, direction_command); 
	vector_sub(3,desired_rate, EKF_bodyrates, rate_command);
//	vector_mult_scalar(rate_command, rate_gain);
	//multiplay rate command by a gain
	//normalize direction_command and multiply by gain 
	//sum the previous 2 steps = output of controller. 

	
}

void sun_transformation(){
	float Ss[3]={0,0,-1}; float Sx=Ss[0]; float Sy=Ss[1];
	float height=0.01; 
	float sqrt=1-powf(Sx,2)-powf(Sy,2);
	float xp=height*Sx/powf(sqrt, 0.5);
	float yp=height*Sy/powf(sqrt, 0.5);
}

void RK4(){
	//Evaluating k1 values 
    Eulers(&I1, &I2, &I3, omega, T, omega_dot);
    float k1_x=omega_dot[0];
	float k1_y=omega_dot[1]; 
	float k1_z=omega_dot[2]; 
	printf("omegax_dot is \t omegay_dot is \t omegaz_dot is \n ");
	printf("%E \t %E \t  %E \n", omega_dot[0], omega_dot[1], omega_dot[2]);
	
	//Evaluating k2 values 
	float state_2[3]; 
	state_2[0]=omega[0]+h*k1_x/2;
	state_2[1]=omega[1]+h*k1_y/2;
	state_2[2]=omega[2]+h*k1_z/2;
	printf("state_2_x is %E \n", state_2[0]);
	printf("state_2_y is %E \n", state_2[1]);
	printf("state_2_z is %E \n", state_2[2]);
	
	Eulers(&I1, &I2, &I3, state_2, T, omega_dot);
    float k2_x=omega_dot[0];
	float k2_y=omega_dot[1]; 
	float k2_z=omega_dot[2]; 
	printf("k2x is \t k2y is \t k2z is \n ");
	printf("%E \t %E \t  %E \n", k2_x, k2_y, k2_z);
	
	//Evaluating k3 values 
	float state_3[3]; 
	state_3[0]=omega[0]+h*k2_x/2;
	state_3[1]=omega[1]+h*k2_y/2;
	state_3[2]=omega[2]+h*k2_z/2;
	printf("state_2_x is %E \n", state_3[0]);
	printf("state_2_y is %E \n", state_3[1]);
	printf("state_2_z is %E \n", state_3[2]);
	
	Eulers(&I1, &I2, &I3, state_3, T, omega_dot);
    float k3_x=omega_dot[0];
	float k3_y=omega_dot[1]; 
	float k3_z=omega_dot[2]; 
    printf("k3x is \t k3y is \t k3z is \n ");
	printf("%E \t %E \t  %E \n", k3_x, k3_y, k3_z);
	
	//Evaluating k4 values 
	float state_4[3];
	state_4[0]=omega[0]+h*k3_x;
    state_4[1]=omega[1]+h*k3_y;
    state_4[2]=omega[2]+h*k3_z;
    Eulers(&I1, &I2, &I3, state_4, T, omega_dot);
    float k4_x=omega_dot[0];
    float k4_y=omega_dot[1];
    float k4_z=omega_dot[2];

    printf("k4x is \t k4y is \t k4z is \n ");
	printf("%E \t %E \t  %E \n", k4_x, k4_y, k4_z);
	
	//Adding it all together 
	float test1=omega[0]+(k1_x+k2_x+k3_x+k4_x)*h/6;
	float test2=omega[1]+(k1_y+k2_y+k3_y+k4_y)*h/6;
	float test3=omega[2]+(k1_z+k2_z+k3_z+k4_z)*h/6;
	printf("omega_x is now %E \n", test1);
	printf("omega_y is now %E \n", test2);
	printf("omega_z is now %E \n", test3);
	
	omega[0]=test1; 
}

//***************************Matrix Math Functions***********************// 
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




//Random print statments I might need 
/*
printf("%E \t %E \t %E \t %E \t %E \n  ", C[0][0], C[0][1], C[0][2], C[0][3],C[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", C[1][0], C[1][1], C[1][2], C[1][3],C[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", C[2][0], C[2][1], C[2][2], C[2][3],C[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", C[3][0], C[3][1], C[3][2], C[3][3],C[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", C[4][0], C[4][1], C[4][2], C[4][3],C[4][4]);/*

printf("Ad is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Jacobian_A[0][0], Jacobian_A[0][1], Jacobian_A[0][2], Jacobian_A[0][3],Jacobian_A[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Jacobian_A[1][0], Jacobian_A[1][1], Jacobian_A[1][2], Jacobian_A[1][3],Jacobian_A[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Jacobian_A[2][0], Jacobian_A[2][1], Jacobian_A[2][2], Jacobian_A[2][3],Jacobian_A[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Jacobian_A[3][0], Jacobian_A[3][1], Jacobian_A[3][2], Jacobian_A[3][3],Jacobian_A[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Jacobian_A[4][0], Jacobian_A[4][1], Jacobian_A[4][2], Jacobian_A[4][3],Jacobian_A[4][4]);

printf("%E \t %E \t %E \n", Qd[0][0], Qd[0][1], Qd[0][2]);
printf("%E \t %E \t %E \n", Qd[1][0], Qd[1][1], Qd[1][2]);
printf("%E \t %E \t %E \n",Qd[2][0], Qd[2][1], Qd[2][2]);
printf("%E \t %E \t %E \n", Qd[3][0], Qd[3][1], Qd[3][2]);
printf("%E \t %E \t %E \n", Qd[4][0], Qd[4][1], Qd[4][2]);

printf("Pk is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[0][0], Pk[0][1], Pk[0][2], Pk[0][3],Pk[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[1][0], Pk[1][1], Pk[1][2], Pk[1][3],Pk[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[2][0], Pk[2][1], Pk[2][2], Pk[2][3],Pk[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[3][0], Pk[3][1], Pk[3][2], Pk[3][3],Pk[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[4][0], Pk[4][1], Pk[4][2], Pk[4][3],Pk[4][4]);

printf("Pk_add_R is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_add_R[0][0], Pk_add_R[0][1], Pk_add_R[0][2], Pk_add_R[0][3],Pk_add_R[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_add_R[1][0], Pk_add_R[1][1], Pk_add_R[1][2], Pk_add_R[1][3],Pk_add_R[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_add_R[2][0], Pk_add_R[2][1], Pk_add_R[2][2], Pk_add_R[2][3],Pk_add_R[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_add_R[3][0], Pk_add_R[3][1], Pk_add_R[3][2], Pk_add_R[3][3],Pk_add_R[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_add_R[4][0], Pk_add_R[4][1], Pk_add_R[4][2], Pk_add_R[4][3],Pk_add_R[4][4]);
*/

/*

	//Evaluating k2 values 
	float state_2[3]; 
	state_2[0]=omega[0]+h*k1_x/2;
	state_2[1]=omega[1]+h*k1_y/2;
	state_2[2]=omega[2]+h*k1_z/2;
	printf("state_2_x is %E \n", state_2[0]);
	printf("state_2_y is %E \n", state_2[1]);
	printf("state_2_z is %E \n", state_2[2]);
	
	Eulers(&I1, &I2, &I3, state_2, T, omega_dot);
    float k2_x=omega_dot[0];
	float k2_y=omega_dot[1]; 
	float k2_z=omega_dot[2]; 
	printf("k2x is \t k2y is \t k2z is \n ");
	printf("%E \t %E \t  %E \n", k2_x, k2_y, k2_z);
	
	//Evaluating k3 values 
	float state_3[3]; 
	state_3[0]=omega[0]+h*k2_x/2;
	state_3[1]=omega[1]+h*k2_y/2;
	state_3[2]=omega[2]+h*k2_z/2;
	printf("state_2_x is %E \n", state_3[0]);
	printf("state_2_y is %E \n", state_3[1]);
	printf("state_2_z is %E \n", state_3[2]);
	
	Eulers(&I1, &I2, &I3, state_3, T, omega_dot);
    float k3_x=omega_dot[0];
	float k3_y=omega_dot[1]; 
	float k3_z=omega_dot[2]; 
    printf("k3x is \t k3y is \t k3z is \n ");
	printf("%E \t %E \t  %E \n", k3_x, k3_y, k3_z);
	
	//Evaluating k4 values 
	float state_4[3];
	state_4[0]=omega[0]+h*k3_x;
    state_4[1]=omega[1]+h*k3_y;
    state_4[2]=omega[2]+h*k3_z;
    Eulers(&I1, &I2, &I3, state_4, T, omega_dot);
    float k4_x=omega_dot[0];
    float k4_y=omega_dot[1];
    float k4_z=omega_dot[2];

    printf("k4x is \t k4y is \t k4z is \n ");
	printf("%E \t %E \t  %E \n", k4_x, k4_y, k4_z);
	
	//Adding it all together 
	float test1=omega[0]+(k1_x+2*k2_x+2*k3_x+k4_x)*h/6;
	float test2=omega[1]+(k1_y+2*k2_y+2*k3_y+k4_y)*h/6;
	float test3=omega[2]+(k1_z+2*k2_z+2*k3_z+k4_z)*h/6;
	printf("omega_x is now %.8E \n", test1);
	printf("omega_y is now %.8E \n", test2);
	printf("omega_z is now %.8E \n", test3);
	
 */      
