#include <stdio.h>
#include <math.h>



void Eulers(float *I1, float *I2, float *I3, float *omega[], float *T[], float *omegax_dot, float *omegay_dot, float *omegaz_dot);


/*variables 
%The value of T is the torque we apply on the spacecraft-this will be calculated, 
for now just leave as a constant*/
//For Euler's Equations 
float T[3] ={1,5,7};
float omega[3]={0.01745, 1, 2}; 
float I1=1; float I2=2; float I3=3; 
float omegax_dot, omegay_dot, omegaz_dot; 
// All omega_previous are initialized by sampling the gyros to read the body rates.
float omegax, omegay, omegaz, omegax_previous, omegay_previous, omegaz_previous;


float h=0.1; // Time step 


int main(){
	

	
	Eulers(I1, I2, I3, omega, T, omegax_dot, omegay_dot, omegaz_dot);
    printf("dot product is %f %f %f", omegax_dot, omegay_dot, omegaz_dot);
	return (0);
}



void Eulers(float *I1, float *I2, float *I3, float *omega[], float *T[], float *omegax_dot, float *omegay_dot, float *omegaz_dot);
    *omegax_dot = I1/(T[0]-(I3-I2)*omega[1]*omega[2]);
	*omegay_dot = I2/(T[1]-(I1-I3)*omega[0]*omega[2]);
	*omegaz_dot = I3/(T[2]-(I2-I1)*omega[0]*omega[1]); 
	
	printf("%f\n", *omegax_dot);
	printf("%f\n", *omegay_dot);
	printf("%f\n", *omegaz_dot);
	
}

