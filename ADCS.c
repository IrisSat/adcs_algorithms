#include<stdio.h>
#include<math.h>
#include "C:\Users\toidjana\Documents\GitHub\adcs_algorithms\matrix_math\matrix_math.h"


float Eulers(float I1, float I2, float I3, float omega[], float T[], float* omegax_dot, float* omegay_dot, float* omegaz_dot);
float RK4(float k1, float k2, float k3, float k4, float h, float previous_estimate);
float quaternions(float omega[], float epsilon[]); 
float Rotation_Matrix(float q[]); 
void dot_prod(float *a, float *b, float *c);


/*variables 
%The value of T is the torque we apply on the spacecraft-this will be calculated, 
for now just leave as a constant*/
float T[3] ={1,5,7};
float omega[3]={0.01745, 1, 2}; 
float I1=1; float I2=2; float I3=3; 
float omegax_dot, omegay_dot, omegaz_dot; 
// All omega_previous are initialized by sampling the gyros to read the body rates.
float omegax, omegay, omegaz, omegax_previous, omegay_previous, omegaz_previous;

//hey 


int main(){
	
	// Step 1: Initialize the omegax_previous, omegy_previous and omegaz_previous to the measurments obtained from the gyro 
	// step 2: Run Euler's function
	// step 3: use multiple outputs of Euler's function to integrate for omegax, omegay, omegaz using RK4 method. 
	// Step 4: Pass the outputs of the integrator to the quaternions function 
	// Step 5: Pass the quaternions to the Rotation_Matrix function 
	// step 6: Initialize the S_ECI vector in body-(Just a random guess)
	// step 7: Use the roation matrix to calculate S_b and S_s from S_ECI
	// step 8: Convert S_s to plane coordinate xp and yp
	// Step 9: Read info from Solar panels to determine if eclipse 
	// Step 10: EKF 
	// Step 11: Convert back to sensor and then body sun measurments
	// step 12: Pass states to controller and calculate output torque 
	// step 13: mapping and saturation 
	// Torque actuation 
	// step 14: Calculate T to feedback into Euler's method 
	// step 15: backup function to "nudge" the spacecraft if it doesn't see the sun for some n amout of time
	Eulers(I1, I2, I3, omega, T, &omegax_dot, &omegay_dot, &omegaz_dot);



	float epsilon[3]={1,2,3};
		
    c_dot=quaternions(omega, epsilon);
    printf("omegax_dot is %f", c_dot);
	return (0);
}

/*
float RK4(float result){
	return(result);
}*/

float Eulers(float I1, float I2, float I3, float omega[], float T[], float* omegax_dot, float* omegay_dot, float* omegaz_dot){
    *omegax_dot = I1/(T[0]-(I3-I2)*omega[1]*omega[2]);
	*omegay_dot = I2/(T[1]-(I1-I3)*omega[0]*omega[2]);
	*omegaz_dot = I3/(T[2]-(I2-I1)*omega[0]*omega[1]); 
	return(*omegax_dot, *omegay_dot, *omegaz_dot);
}



float RK4(float k1, float k2, float k3, float k4, float h, float previous_estimate){
	// Numerical integration using the Runge Kutta 4 method 
	// To calculate the rk4 method, we need to sample we need 3 data points from Euler's equation, seperated by some time step h/2 
	//k1=f(t), k2=k3=f(t+h/2), k4=f(t+h)
	float apriori_estimate=previous_estimate+(1/6)*h*(k1+k2+k3+k4);
	return apriori_estimate; 
}


float quaternions(float omega[], float epsilon[]){
	// m rows, n coloumns 
	float c_dot; 
	dot_prod(omega,epsilon, &c_dot);
	return &c_dot; 
	
}
/*
float Rotation_Matrix(q[]){
	
}*/



//Matrix Math functions 
void dot_prod(float *a, float *b, float *c)
{
    *c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
