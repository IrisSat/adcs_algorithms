// ADCS algorithms Version 4
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void Calc_Jacobian(float *I1, float *I2, float *I3, float *height, float *omega, float *S_p, float Jacobian_A[][5]);
void CalcQd(float *I1, float *I2, float *I3, float *w, float Qd[][5]);
void PredictP(float Jacobian_A[][5], float Pk1[][5], float Qd[][5]);
void CalcK(float Pk1[][5], float R[][5], float K[][5]); 
void RK4(float I1, float I2, float I3, float *omega, float *T);
void EKF(float I1, float I2, float I3, float *omega, float *EKF_state);
void calc_omegax_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegax_dot );
void calc_omegay_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegay_dot );
void calc_omegaz_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegaz_dot );
void DipoleMapping(float *commad_torque, float *B_field, float *command_dipole);
void saturation(float *command_dipole);
void sun_transform(float *EKF_planar, int eclipse_flag, int front_panel_flag, float *S_body); 
void controller(float *bodyrates, float *S_body, float *torque_command);
void reorient(float *B_field, float *EKF_bodyrates);  

//***************************Matrix Math Functions***********************//  
//Note I changed some of the functions below from their original version in matrix_math file. 
void matmul(int m, int n, int p, float A[m][n], float B[n][p], float C[m][p]);
void transpose(int m, int n, float A[][n], float A_T[][m]);
void addmatrix(int m, int n, float A[][n], float B[][n], float C[][n]);
void vectorMultiply(int m, int n, float A[m][n], float B[m], float C[m]);
void vector_sub(int m, float *a, float *b, float *c);
void cross_prod(float *a, float *b, float *c);
void vector_add(int m, float *a, float *b, float *c);
void invert(int m, float A[m][m], float AInv[m][m]);
void solveInv(int m, float A[m][m], float B[m][m], float X[m][m]);
void dot_prod(float *a, float *b, float *c);
void vector_mult_scalar(float *v, float mult);
void vector_div_scalar(float *v, float div); 
void norm(float *a, float *normA);
//************************************************************************//

// Variables required for Euler's Equations 
float I1=3.6287e-2; float I2=3.4092e-2; float I3=0.9026e-2; 
float omega_dot[3];
float omegax_dot; float omegay_dot; float omegaz_dot; float torque_command[3]; 
float T[3] ={0,0,0};


int eclipse_flag=0; 
//when eclipse flag is 1, the spacecraft is in eclipse. 
int front_panel_flag=1; 
//when front_panel_flag is 1, the front solar panels are facing the sun-should have a current reading. 


int main(){
	
//Step 1: sample Gyros to obtain a initial rate for integration 
float omega[3]={0.01745, 0.001745, 0.001745}; 
float EKF_state[5]; float EKF_bodyrates[3]; float EKF_planar[2]; 

//From here, start a loop 
//printf("omega is [%.8E %.8E %.8E]\n ", omega[0], omega[1], omega[2] ); 
RK4(I1, I2, I3, omega, T);
//printf("omega is [%.8E %.8E %.8E]\n ", omega[0], omega[1], omega[2] ); 
EKF(I1, I2, I3, omega, EKF_state);	
EKF_bodyrates[0]=EKF_state[0]; EKF_bodyrates[1]=EKF_state[1]; EKF_bodyrates[2]=EKF_state[2]; 
EKF_planar[0]=EKF_state[3]; EKF_planar[1]=EKF_state[4]; 

//Convert planar measurments to sun vector 
float planar[2]={-.1, 0.6}; 
float S_body[3]; 
sun_transform(planar, eclipse_flag, front_panel_flag, S_body);
//printf("S_body after solar panel check is %f %f %f \n", S_body[0], S_body[1],S_body[2]); 

//Controller 
controller(EKF_bodyrates, S_body, torque_command); 

//Dipole mapping 
float command_dipole[3];
float B_field[3]={0.1892, 0.005, -0.013};  //This will be a magnetometer measurment 
DipoleMapping(torque_command, B_field, command_dipole); 

//Saturation 
saturation(command_dipole); 

//Calculate appplied torque and set equal to T 
cross_prod(command_dipole, B_field, T); 
//printf("T is %.8E %.8E %.8E \n", T[0], T[1],T[2]); 

}

void controller(float *bodyrates, float *S_body, float *torque_command){
	//The controller takes in the estimated body rates and sun measurments determined by the EKF.
	//The controller outputs a command torque
	float desired_direction[3]={-1,0,0}; 
		printf("body rates are is %.8E %.8E %.8E \n", bodyrates[0], bodyrates[1], bodyrates[2] );

	float desired_rate[3]={0.01745329, 0, 0}; 
	float rate_gain=0.0004; float pointing_gain=-0.000005; 
	float pointing_command[3]; float rate_command[3]; float pointing_cross_prod[3];
	cross_prod(S_body, desired_direction, pointing_cross_prod); 
	vector_sub(3,desired_rate, bodyrates, rate_command);
	printf("body rates are is %.8E %.8E %.8E \n", bodyrates[0], bodyrates[1], bodyrates[2] );

	printf("rate command is %.8E %.8E %.8E \n", rate_command[0], rate_command[1], rate_command[2] );

	vector_mult_scalar(rate_command, rate_gain);          
	norm(pointing_cross_prod, pointing_command); 
	vector_mult_scalar(pointing_command, pointing_gain);
	vector_add(3, rate_command, pointing_command, torque_command);
	printf("rate command is %.8E %.8E %.8E \n", rate_command[0], rate_command[1], rate_command[2] );
	printf("Torque command is %.8E %.8E %.8E \n", torque_command[0], torque_command[1], torque_command[2] );
}

void DipoleMapping(float *commad_torque, float *B_field, float *command_dipole){
	float B_dot=0; 
	dot_prod(B_field, B_field, &B_dot); 
	cross_prod(B_field, commad_torque, command_dipole); 
	vector_div_scalar(command_dipole, B_dot);
	//printf("m is %.8E %.8E %.8E \n", command_dipole[0], command_dipole[1], command_dipole[2]);
}


void sun_transform(float *EKF_planar, int eclipse_flag, int front_panel_flag, float *S_body){
	float xp=EKF_planar[0]; 
	float yp=EKF_planar[1]; 
	float height=0.01; 
	float sx=sqrt(1/(1+(yp*yp)/(xp*xp)+(height*height)/(xp*xp)));  //sx=sqrt(1/(1+yp^2/xp^2+h^2/xp^2));
    if (xp <0){
    	sx=-sx; 
	}
	float sy=yp/xp*sx; 
	float sz=sqrt(1-(sx*sx)-(sy*sy));
	//printf("Sx sy sz is %f %f %f \n", sx, sy,sz); 

//This is equivalent to a 90 degree 2nd axis rotation. 
	S_body[0]=sz;
	S_body[1]=sy;
	S_body[2]=-sx; 
	//printf("S_body is %f %f %f \n", S_body[0], S_body[1],S_body[2]); 
	
//This is using the solar cell input 
	if (front_panel_flag ==1){   //if there is sun on the front solar panel 
		S_body[0]=-S_body[0]; 
	}
	if (eclipse_flag == 1){      //if the spacecraft is in eclipse 
	    S_body[0]=-S_body[0];
    }
	

}


void saturation(float *command_dipole){
	float limit=0.1731;
//Take absolute value of command_dipole array 
    float command_dipole_new[3]={fabs(command_dipole[0]), fabs(command_dipole[1]), fabs(command_dipole[2])} ;

//Finding the maximum value
    int i; float max;
    for (i=0; i<3; i++){
    	if (max < command_dipole_new[i]){
    	    max=command_dipole_new[i];
    	
    	}
	}	
//	printf("The maximum in the array is %.8E \n", max);
	
	if (max>=limit){
		//If one of the elements in the command dipole array is larger than limit, we scale the whole command dipole vector 
		//If all elements are smaller, command dipole stays as is. 
	    vector_mult_scalar(command_dipole, (limit/max));
	}


}

void EKF(float I1, float I2, float I3, float *omega, float *EKF_state){
//Will need to make this take in the true state propogation. 
// Variables required for EKF 
float w=1e-5; //Assumed process noise 
float height=0.01; 

float Pk1[5][5]={
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
PredictP(Jacobian_A,  Pk1, Qd);
CalcK(Pk1, R, K);

float True_state[5]={omega[0],omega[1],omega[2],S_p[0],S_p[1]} ; //Move this 
float Measured_state[5]={0.0221420309324296, 0.0177486658038343, -0.0179671575114262, 0.369903257448045, -0.380396852672751};

// Check for eclipse 
if (eclipse_flag ==1){
	H[3][3]=0;
	H[4][4]=0; 
}

//State updates 
float H_times_True_state[5];
vectorMultiply(5, 5, H, True_state, H_times_True_state);
float M_minus_T[5];
float K_times_M_minus_T[5]; 
vector_sub(5, Measured_state, H_times_True_state, M_minus_T);
vectorMultiply(5,5, K, M_minus_T, K_times_M_minus_T );
vector_add(5, True_state, K_times_M_minus_T, EKF_state  );




//printf("%E \t %E \t %E \t %E \t %E \n  ", Xk[0], Xk[1], Xk[2], Xk[3],Xk[4]);
	
}



void calc_omegax_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegax_dot ){
	*omegax_dot=(T[0]-(*I3-*I2)*omega[1]*omega[2])/(*I1);
}

void calc_omegay_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegay_dot ){
	*omegay_dot=(T[1]-(*I1-*I3)*omega[0]*omega[2])/(*I2);
}

void calc_omegaz_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegaz_dot ){
	*omegaz_dot=(T[2]-(*I2-*I1)*omega[0]*omega[1])/(*I3); 
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


void PredictP(float Jacobian_A[][5], float Pk1[][5], float Qd[][5]){
	// This function calculates & updates the error covariance (Pk)
	float C[5][5];
	float D[5][5];
    float E[5][5];
    
	matmul(5,5,5,Jacobian_A, Pk1, C);
	transpose(5,5, Jacobian_A,D);
	matmul(5,5,5,C,D,E);
	addmatrix(5,5,E,Qd,Pk1);

}

void CalcK(float Pk1[][5], float R[][5], float K[][5]){
	
//The equation to calculate K simplifies to Pk/(Pk+R) if H is just the identity matrix. 
float Pk_add_R[5][5]; 
addmatrix(5, 5, Pk1, R, Pk_add_R );
solveInv(5, Pk_add_R, Pk1, K);
/*
printf("K is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", K[0][0], K[0][1], K[0][2], K[0][3],K[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[1][0], K[1][1], K[1][2], K[1][3],K[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[2][0], K[2][1], K[2][2], K[2][3],K[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[3][0], K[3][1], K[3][2], K[3][3],K[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", K[4][0], K[4][1], K[4][2], K[4][3],K[4][4]);	*/
}

void RK4(float I1, float I2, float I3, float *omega, float *T){
    float h=0.1; //Step size (s)
    float omega_dot; float omegax; float omegay; float omegaz; 
    float k1_x; float k2_x; float k3_x; float k4_x;
    float k1_y; float k2_y; float k3_y; float k4_y;
    float k1_z; float k2_z; float k3_z; float k4_z;
//Integrating omega x 
// k1 value 
	calc_omegax_dot(&I1, &I2, &I3, omega, T, &omega_dot);  
    k1_x=omega_dot;
   
// K2 value 
    float omega_new[3]={omega[0]+k1_x*h/2, omega[1], omega[2]};
    calc_omegax_dot(&I1, &I2, &I3, omega_new, T, &omega_dot);
    k2_x=omega_dot; 

// k3 value 
    omega_new[0]=omega[0]+k2_x*h/2;
    calc_omegax_dot(&I1, &I2, &I3, omega_new, T, &omega_dot);
    k3_x=omega_dot; 

// k4 value 
    omega_new[0]=omega[0]+k3_x*h;
    calc_omegax_dot(&I1, &I2, &I3, omega_new, T, &omega_dot);
    k4_x=omega_dot; 

//Putting it all together for omega_x
	omegax=omega[0]+(k1_x+2*k2_x+2*k3_x+k4_x)*h/6;
	
//__________________________________________________________________________________//	
//Integrating omega y 
// k1 value 
	calc_omegay_dot(&I1, &I2, &I3, omega, T, &omega_dot);  
    k1_y=omega_dot;
   
// K2 value 
    omega_new[0]=omega[0]; omega_new[1]=omega[1]+k1_y*h/2;  omega_new[2]=omega[2]; 
    calc_omegay_dot(&I1, &I2, &I3, omega_new, T, &omega_dot);
    k2_y=omega_dot; 

// k3 value 
    omega_new[1]=omega[1]+k2_y*h/2;
    calc_omegay_dot(&I1, &I2, &I3, omega_new, T, &omega_dot);
    k3_y=omega_dot; 

// k4 value 
    omega_new[1]=omega[1]+k3_y*h;
    calc_omegay_dot(&I1, &I2, &I3, omega_new, T, &omega_dot);
    k4_y=omega_dot; 

//Putting it all together for omega_x
	omegay=omega[1]+(k1_y+2*k2_y+2*k3_y+k4_y)*h/6;
		
//__________________________________________________________________________________//	
//Integrating omega z 
// k1 value 
	calc_omegaz_dot(&I1, &I2, &I3, omega, T, &omega_dot);  
    k1_z=omega_dot;
   
// K2 value 
    omega_new[0]=omega[0]; omega_new[1]=omega[1];  omega_new[2]=omega[2]+k1_z*h/2; 
    calc_omegaz_dot(&I1, &I2, &I3, omega_new, T, &omega_dot);
    k2_z=omega_dot; 

// k3 value 
    omega_new[2]=omega[2]+k2_z*h/2;
    calc_omegaz_dot(&I1, &I2, &I3, omega_new, T, &omega_dot);
    k3_z=omega_dot; 

// k4 value 
    omega_new[2]=omega[2]+k3_z*h;
    calc_omegaz_dot(&I1, &I2, &I3, omega_new, T, &omega_dot);
    k4_z=omega_dot; 

//Putting it all together for omega_x
	omegaz=omega[2]+(k1_z+2*k2_z+2*k3_z+k4_z)*h/6;
	
//Rewrite omega to new calculated values 
omega[0]=omegax; 
omega[1]=omegay; 
omega[2]=omegaz; 	
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

void vector_add(int m, float *a, float *b, float *c)
{
    for(int i=0; i < m; i++)
        c[i] = a[i]+b[i];
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

void dot_prod(float *a, float *b, float *c)
{
    *c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void vector_div_scalar(float *v, float div)
{
    for(int i=0; i < 3; i++)
        v[i] /= div;
}

void vector_mult_scalar(float *v, float mult)
{
    for(int i=0; i < 3; i++)
        v[i] *= mult;
}

void norm(float *a, float *normA){
    double mag = sqrt( (a[0]*a[0]) + (a[1]*a[1]) + (a[2]*a[2]));
    normA[0] = a[0] / mag;
    normA[1] = a[1] / mag;
    normA[2] = a[2] / mag;
}