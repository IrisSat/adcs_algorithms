#include <stdio.h>
#include <math.h>

void Eulers(float *I1, float *I2, float *I3, float *omega, float *T, float *omega_dot );
void RK4();



float I1=3.6287e-2; float I2=3.4092e-2; float I3=0.9026e-2; 
float omega_dot[3];
float T[3] ={0,0,0};
float omega[3]={0.01745, 0.001745, 0.001745}; 
float h=0.1; //Step size (s)







int main(){
RK4();
printf("%E ", omega[0]);
}


void Eulers(float *I1, float *I2, float *I3, float *omega, float *T, float *omega_dot ){
	omega_dot[0]=(T[0]-(*I3-*I2)*omega[1]*omega[2])/(*I1);
	omega_dot[1]=(T[1]-(*I1-*I3)*omega[0]*omega[2])/(*I2);
    omega_dot[2]=(T[2]-(*I2-*I1)*omega[0]*omega[1])/(*I3); 
//	printf("omegax_dot is \t omegay_dot is \t omegaz_dot is \n ");
//	printf("%E \t %E \t  %E", omega_dot[0], omega_dot[1], omega_dot[2]);
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