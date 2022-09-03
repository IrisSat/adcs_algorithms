#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void Calc_Jacobian(float *I1, float *I2, float *I3, float *height, float *omega, float *S_p, float Jacobian_A[][5]);
void CalcQd(float *I1, float *I2, float *I3, float *w, float Qd[][5]);
void PredictP(float Jacobian_A[][5], float Pk[][5], float Qd[][5]);
void CalcK(float Pk[][5], float R[][5], float K[][5], float H[][5]); 
void CorrectP(float Pk_old[][5], float K[][5], float H[][5], float Pk_corrected[][5]);
void RK4(float I1, float I2, float I3, float *omega, float *T);
void EKF(float I1, float I2, float I3, float Pk[][5], float *omega, float *EKF_state, float *Gyro_read, float *Sun_read, int eclipse_flag);
void calc_omegax_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegax_dot );
void calc_omegay_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegay_dot );
void calc_omegaz_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegaz_dot );
void DipoleMapping(float *commad_torque, float *B_field, float *command_dipole);
void saturation(float *command_dipole);
void sun_transform(float *EKF_planar, int eclipse_flag, int front_panel_flag, float *S_body); 
void controller(float *bodyrates, float *S_body, float *torque_command);
void reorient(float *B_field, float *command_dipole);

//***************************Matrix Math Functions***********************//  
//Note I changed some of the functions below from their original version in matrix_math file. 
void matmul(int m, int n, int p, float A[m][n], float B[n][p], float C[m][p]);
void transpose(int m, int n, float A[][n], float A_T[][m]);
void addmatrix(int m, int n, float A[][n], float B[][n], float C[][n]);
void submatrix(int m, int n, float A[][n], float B[][n], float C[][n]);

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




int main(){
	
	//Initialize variables 
	float Gyro_read[3]; 
	float Mag_read[3]; 
	float Sun_read[2];
	int eclipse_flag;
	int front_panel_flag;
	
    //During the first iteration, initialize omega to the gyro measurments. 
	float omega[3]={0.01745, 0.001745, 0.001745}; //These are random numbers but should be gyro readings
	float EKF_state[5]; 
	//Initialize EKF_state. The first 3 values get initiated to the same gyro measurments that were used to initiate omega
	EKF_state[0]=omega[0]; EKF_state[1]=omega[1]; EKF_state[2]=omega[2]; EKF_state[3]=0; EKF_state[4]=0; 
	float EKF_bodyrates[3]; float EKF_planar[2]; 
	
	float I1=3.6287e-2; float I2=3.4092e-2; float I3=0.9026e-2; 
	
	//Initial value of Pk
   	float Pk[5][5]={
     {0.0076154,0,0,0,0}, 
     {0,0.0076154,0,0,0}, 
     {0,0,0.0076154,0,0}, 
     {0,0,0,304.617,0}, 
     {0,0,0,0,304.617}
     };
      
    float S_body[3]; //Output variable of Sun_transofrm function 
    float torque_command[3]; 
    float command_dipole[3];
    float T[3] ={0,0,0};


	
for  (;;){
	
    //Sample gyros and sun sensors, obtain information about the solar cell currents.  
     Gyro_read[0]=0.01745;  Gyro_read[1]=0.001745;  Gyro_read[2]=0.001745;     //Sample Gyros 
     Sun_read[0]=-.1; Sun_read[1]=0.6;                                         //Sample Sunsensors
	 eclipse_flag=0;                                                           //Check eclipse flag, when eclipse flag is 1, the spacecraft is in eclipse. 
     front_panel_flag=1;                                                       //Check front panels, when front_panel_flag is 1, the front solar panels are facing the sun-should have a current reading. 
     
	  
      RK4(I1, I2, I3, omega, T);  //Propogates omega ahead of time. Re-writes the omega variable. 
      
      EKF(I1, I2, I3, Pk, omega, EKF_state, Gyro_read, Sun_read, eclipse_flag);	
      EKF_bodyrates[0]=EKF_state[0]; EKF_bodyrates[1]=EKF_state[1]; EKF_bodyrates[2]=EKF_state[2]; 
      EKF_planar[0]=EKF_state[3]; EKF_planar[1]=EKF_state[4]; 

     //Convert sun sensor measurments from EKF output to sun vector 
     sun_transform(EKF_planar, eclipse_flag, front_panel_flag, S_body);
    //printf("S_body after solar panel check is %f %f %f \n", S_body[0], S_body[1],S_body[2]); 

    //Controller 
     controller(EKF_bodyrates, S_body, torque_command);                       //Calculates the torque command based on the estimated body rates and sun sensor position. 

    //Dipole mapping 
    Mag_read[0]=0.001; Mag_read[1]=-0.0023; Mag_read[2]=-0.007;               //Sample the magnetometers 
    DipoleMapping(torque_command, Mag_read, command_dipole);                  //Maps the required torque command to a dipole 

    //Saturation 
    saturation(command_dipole);                                               //Scales the command dipole. The command_dipole variable gets re-written. 

    //Calculate appplied torque and set equal to T                            //This T gets used in the RK4 in the next iteration. 
    cross_prod(command_dipole, Mag_read, T); 
    //printf("T is %.8E %.8E %.8E \n", T[0], T[1],T[2]); 

    //Actuate torque rods with command dipole 

}
//reorient\plan B 

//reorient(B_field, command_dipole); 
//printf("command dipole from reorient function is %.8E %.8E %.8E \n", command_dipole[0], command_dipole[1],command_dipole[2]); 

}
void reorient(float *B_field, float *command_dipole){
	//Activate this function if the sun sensors have not been able to detect the sun for more than 1 orbit. 
	float limit=0.1731; 
    float By=B_field[1]; float Bz=B_field[2]; 
	
	command_dipole[0]=limit/2; 
	command_dipole[1]=limit/2; 
	command_dipole[2]=command_dipole[1]*Bz/By; 
	//Hold the command dipole for 10 seconds. 
	//Activate again if sun is not detected by the sun sensors for more than 1 orbit. 	
	//Make sure the spacecraft returns to its normal control operations as soon as this is done. 
	//Nudging the spacecraft this way will cause an increase in the body rates that we will want to control as soon as possible to prevent an unstable spin. 
}




void controller(float *bodyrates, float *S_body, float *torque_command){
	//The controller takes in the estimated body rates and sun measurments determined by the EKF.
	//The controller outputs a command torque
	float desired_direction[3]={-1,0,0}; 
	float desired_rate[3]={0.01745329, 0, 0}; 
	float rate_gain=0.0004; float pointing_gain=-0.000005; 
	float pointing_command[3]; float rate_command[3]; float pointing_cross_prod[3];
	
	cross_prod(S_body, desired_direction, pointing_cross_prod); 
	vector_sub(3,desired_rate, bodyrates, rate_command);
	vector_mult_scalar(rate_command, rate_gain);          
	norm(pointing_cross_prod, pointing_command); 
	vector_mult_scalar(pointing_command, pointing_gain);
	vector_add(3, rate_command, pointing_command, torque_command);

}

void DipoleMapping(float *commad_torque, float *B_field, float *command_dipole){
	//This function maps the require torque command to a possible dipole based on the magnetic field measurments. 
	float B_dot=0; 
	dot_prod(B_field, B_field, &B_dot); 
	cross_prod(B_field, commad_torque, command_dipole); 
	vector_div_scalar(command_dipole, B_dot);
	//printf("m is %.8E %.8E %.8E \n", command_dipole[0], command_dipole[1], command_dipole[2]);
}


void sun_transform(float *EKF_planar, int eclipse_flag, int front_panel_flag, float *S_body){
	//This function transforms the sun sensor measurment from 2D to a 3D vector. S_body is the 3D vector
	float xp=EKF_planar[0]; 
	float yp=EKF_planar[1]; 
	float height=0.01;  //WILL NEED TO UPDATE THIS BASED ON SUN SENSOR DATA SHEET
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
	//Scales the command dipole based on the operating range of the rods. command_dipole is rewritten in this funcion. 
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
	if (max>=limit){
		//If one of the elements in the command dipole array is larger than limit, we scale the whole command dipole vector 
		//If all elements are smaller, command dipole stays as is. 
	    vector_mult_scalar(command_dipole, (limit/max));
	}
}

void EKF(float I1, float I2, float I3, float Pk[][5], float *omega, float *EKF_state, float *Gyro_read, float *Sun_read, int eclipse_flag){

float w=1e-5;        //Assumed process noise 
float height=0.01;   //WILL NEED TO UPDATE THIS BASED ON SUN SENSOR DATA SHEET
float Qd[5][5];	
float K [5][5];
float S_p[2]={0,0};  

float True_state[5]={omega[0],omega[1],omega[2],S_p[0],S_p[1]} ; 
//float Measured_state[5]={0.0221420309324296, 0.0177486658038343, -0.0179671575114262, 0.369903257448045, -0.380396852672751}; Used for testing.
float Measured_state[5]; 
Measured_state[0]=Gyro_read[0];  Measured_state[1]=Gyro_read[1];  Measured_state[2]=Gyro_read[2];   //Complete this 
Measured_state[3]=Sun_read[0]; Measured_state[4]=Sun_read[1];


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

float Jacobian_A[5][5]={
{0,0,0,0,0},
{0,0,0,0,0},
{0,0,0,0,0},
{0,0,0,0,0},
{0,0,0,0,0}
};


/*
//For testing, delete when done. 
omega[0]=0.01745;
omega[1]= 0.001745;
omega[2]=0.001745;
S_p[0]=0.34185;
S_p[1]=0.35156;
*/

//Calc jacobian is calculated based on the latest EKF estimate 
float body_rates_previous[3]; float sun_previous[2]; 
body_rates_previous[0]=EKF_state[0]; body_rates_previous[1]=EKF_state[1]; body_rates_previous[2]=EKF_state[2];
sun_previous[0]=EKF_state[3]; sun_previous[1]=EKF_state[3];

Calc_Jacobian(&I1, &I2, &I3, &height, body_rates_previous, sun_previous, Jacobian_A);
CalcQd(&I1, &I2, &I3, &w, Qd);
PredictP(Jacobian_A,  Pk, Qd);
/*
printf("PredictedP \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[0][0], Pk[0][1], Pk[0][2], Pk[0][3],Pk[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[1][0], Pk[1][1], Pk[1][2], Pk[1][3],Pk[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[2][0], Pk[2][1], Pk[2][2], Pk[2][3],Pk[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[3][0], Pk[3][1], Pk[3][2], Pk[3][3],Pk[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[4][0], Pk[4][1],Pk[4][2], Pk[4][3],Pk[4][4]);
*/
// Check for eclipse 


if (eclipse_flag ==1){
	H[3][3]=0;
	H[4][4]=0; 
}


CalcK(Pk, R, K, H );
/*
printf("K is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", K[0][0], K[0][1], K[0][2], K[0][3],K[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[1][0], K[1][1], K[1][2], K[1][3],K[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[2][0], K[2][1], K[2][2], K[2][3],K[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[3][0], K[3][1], K[3][2], K[3][3],K[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", K[4][0], K[4][1], K[4][2], K[4][3],K[4][4]);
*/
CorrectP(Pk, K, H, Pk);

/*
printf("Pk_corrected in EKF \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[0][0], Pk[0][1], Pk[0][2], Pk[0][3],Pk[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[1][0], Pk[1][1], Pk[1][2], Pk[1][3],Pk[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[2][0], Pk[2][1], Pk[2][2], Pk[2][3],Pk[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[3][0], Pk[3][1], Pk[3][2], Pk[3][3],Pk[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[4][0], Pk[4][1],Pk[4][2], Pk[4][3],Pk[4][4]);	
*/



//State updates 
float H_times_True_state[5];
vectorMultiply(5, 5, H, True_state, H_times_True_state);
float M_minus_T[5];
float K_times_M_minus_T[5]; 
vector_sub(5, Measured_state, H_times_True_state, M_minus_T);
vectorMultiply(5,5, K, M_minus_T, K_times_M_minus_T );
vector_add(5, True_state, K_times_M_minus_T, EKF_state  );
}





void Calc_Jacobian(float *I1, float *I2, float *I3, float *height, float *omega, float *S_p, float Jacobian_A[][5]){	
//Calculates the Jacobian based on the latest estimate of the state. 
float xp=S_p[0]; float yp=S_p[1]; 
float w1=omega[0]; float w2=omega[1]; float w3=omega[2]; 

Jacobian_A[0][1]=-((*I3-*I2)/(*I1))*w3;
Jacobian_A[0][2]=-((*I3-*I2)/(*I1))*w2;

Jacobian_A[1][0]=-((*I1-*I3)/(*I2))*w3; 
Jacobian_A[1][2]=-((*I1-*I3)/(*I2))*w1; 

Jacobian_A[2][0]=-((*I2-*I1)/(*I3))*w2;
Jacobian_A[2][1]=-((*I2-*I1)/(*I3))*w1;

Jacobian_A[3][0]=-yp*xp/(*height);
Jacobian_A[3][1]=(*height)+(xp*xp)/(*height); 
Jacobian_A[3][2]=-yp;
Jacobian_A[3][3]=(-yp*w1/(*height))+(2*xp*w2/(*height));
Jacobian_A[3][4]= -w3-(xp*w1/(*height));

Jacobian_A[4][0]=-(*height)-(yp*yp/(*height)); 
Jacobian_A[4][1]=xp*yp/(*height); 
Jacobian_A[4][2]=xp;
Jacobian_A[4][3]=w3+(yp*w2/(*height));
Jacobian_A[4][4]=(2*yp*w1)/(*height)+(xp*w2/(*height));

/*
printf("Jacobian_A is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Jacobian_A[0][0], Jacobian_A[0][1], Jacobian_A[0][2], Jacobian_A[0][3],Jacobian_A[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Jacobian_A[1][0], Jacobian_A[1][1], Jacobian_A[1][2], Jacobian_A[1][3],Jacobian_A[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Jacobian_A[2][0], Jacobian_A[2][1], Jacobian_A[2][2], Jacobian_A[2][3],Jacobian_A[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Jacobian_A[3][0], Jacobian_A[3][1], Jacobian_A[3][2], Jacobian_A[3][3],Jacobian_A[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Jacobian_A[4][0], Jacobian_A[4][1], Jacobian_A[4][2],Jacobian_A[4][3],Jacobian_A[4][4]);	
*/

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

}

void CorrectP(float Pk_old[][5], float K[][5], float H[][5], float Pk_corrected[][5]){
float I[5][5]={
{1,0,0,0,0},
{0,1,0,0,0},
{0,0,1,0,0},
{0,0,0,1,0},
{0,0,0,0,1}
};

float A[5][5]; 
float B[5][5]; 

	matmul(5,5,5,K, H, A);
	submatrix(5,5,I, A, B );
/*	
printf("A \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", A[0][0], A[0][1], A[0][2], A[0][3],A[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", A[1][0], A[1][1], A[1][2], A[1][3],A[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", A[2][0], A[2][1], A[2][2], A[2][3],A[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", A[3][0], A[3][1], A[3][2], A[3][3],A[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", A[4][0], A[4][1],A[4][2], A[4][3],A[4][4]);
*/	
	matmul(5,5,5,B, Pk_old,Pk_corrected);
	
/*	
printf("Pk_corrected is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_corrected[0][0], Pk_corrected[0][1], Pk_corrected[0][2], Pk_corrected[0][3],Pk_corrected[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk_corrected[1][0], Pk_corrected[1][1], Pk_corrected[1][2], Pk_corrected[1][3],Pk_corrected[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk_corrected[2][0], Pk_corrected[2][1], Pk_corrected[2][2], Pk_corrected[2][3],Pk_corrected[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk_corrected[3][0], Pk_corrected[3][1], Pk_corrected[3][2], Pk_corrected[3][3],Pk_corrected[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_corrected[4][0], Pk_corrected[4][1],Pk_corrected[4][2], Pk_corrected[4][3],Pk_corrected[4][4]);	*/
}


void CalcK(float Pk[][5], float R[][5], float K[][5], float H[][5]){
	
//The equation to calculate K simplifies to Pk/(Pk+R) if H is just the identity matrix. 
float Pk_add_R[5][5]; 
float H_transpose[5][5]; 
float Pk_times_H_transpose[5][5]; 
float H_times_Pk_times_H_transpose[5][5]; 

transpose(5,5,H,H_transpose);
matmul(5,5,5, Pk, H_transpose, Pk_times_H_transpose);
matmul(5,5,5,H,Pk_times_H_transpose, H_times_Pk_times_H_transpose);
addmatrix(5, 5, H_times_Pk_times_H_transpose, R, Pk_add_R );
/*
 printf("Pk_times_H_transpose \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_times_H_transpose[0][0], Pk_times_H_transpose[0][1], Pk_times_H_transpose[0][2], Pk_times_H_transpose[0][3],Pk_times_H_transpose[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk_times_H_transpose[1][0], Pk_times_H_transpose[1][1], Pk_times_H_transpose[1][2], Pk_times_H_transpose[1][3],Pk_times_H_transpose[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk_times_H_transpose[2][0], Pk_times_H_transpose[2][1], Pk_times_H_transpose[2][2], Pk_times_H_transpose[2][3],Pk_times_H_transpose[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk_times_H_transpose[3][0], Pk_times_H_transpose[3][1], Pk_times_H_transpose[3][2], Pk_times_H_transpose[3][3],Pk_times_H_transpose[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_times_H_transpose[4][0], Pk_times_H_transpose[4][1], Pk_times_H_transpose[4][2], Pk_times_H_transpose[4][3],Pk_times_H_transpose[4][4]);

 printf("Pk_add_R is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_add_R[0][0], Pk_add_R[0][1], Pk_add_R[0][2], Pk_add_R[0][3],Pk_add_R[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk_add_R[1][0], Pk_add_R[1][1], Pk_add_R[1][2], Pk_add_R[1][3],Pk_add_R[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk_add_R[2][0], Pk_add_R[2][1], Pk_add_R[2][2], Pk_add_R[2][3],Pk_add_R[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk_add_R[3][0], Pk_add_R[3][1], Pk_add_R[3][2], Pk_add_R[3][3],Pk_add_R[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk_add_R[4][0], Pk_add_R[4][1], Pk_add_R[4][2], Pk_add_R[4][3],Pk_add_R[4][4]);	
solveInv(5, Pk_add_R, Pk_times_H_transpose, K);

/*
printf("Pk_corrected in EKF \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[0][0], Pk[0][1], Pk[0][2], Pk[0][3],Pk[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[1][0], Pk[1][1], Pk[1][2], Pk[1][3],Pk[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[2][0], Pk[2][1], Pk[2][2], Pk[2][3],Pk[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", Pk[3][0], Pk[3][1], Pk[3][2], Pk[3][3],Pk[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", Pk[4][0], Pk[4][1],Pk[4][2], Pk[4][3],Pk[4][4]);	




printf("K is \n");
printf("%E \t %E \t %E \t %E \t %E \n  ", K[0][0], K[0][1], K[0][2], K[0][3],K[0][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[1][0], K[1][1], K[1][2], K[1][3],K[1][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[2][0], K[2][1], K[2][2], K[2][3],K[2][4]);
printf("%E \t %E \t %E \t %E \t %E \n ", K[3][0], K[3][1], K[3][2], K[3][3],K[3][4]);
printf("%E \t %E \t %E \t %E \t %E \n  ", K[4][0], K[4][1], K[4][2], K[4][3],K[4][4]);	
*/
}

//Calculate Eulers equations in every axis-used in the RK4 function. 
void calc_omegax_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegax_dot ){
	*omegax_dot=(T[0]-(*I3-*I2)*omega[1]*omega[2])/(*I1);
}

void calc_omegay_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegay_dot ){
	*omegay_dot=(T[1]-(*I1-*I3)*omega[0]*omega[2])/(*I2);
}

void calc_omegaz_dot(float *I1, float *I2, float *I3, float *omega, float *T, float *omegaz_dot ){
	*omegaz_dot=(T[2]-(*I2-*I1)*omega[0]*omega[1])/(*I3); 
}

void RK4(float I1, float I2, float I3, float *omega, float *T){
    float h=0.1; //Step size (s) //THIS WILL NEED TO BE UPDATED. 
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

void submatrix(int m, int n, float A[][n], float B[][n], float C[][n])
{
    int i,j;
    for(i=0; i < m; i++){
        for(j=0; j < n; j++){
            C[j][i] = A[j][i]-B[j][i];
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