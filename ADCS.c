#include<stdio.h>
#include<math.h>


float Eulers(float I1, float I2, float I3, float omega[], float T[], float* omegax_dot, float* omegay_dot, float* omegaz_dot);
float RK4(float k1, float k2, float k3, float k4, float h, float previous_estimate);
void quaternions(float omega[], float epsilon[], float eta, float* eta_dot, float *epsilon_dot); 
void Rotation_Matrix_Body(float epsilon[], float eta, float C_b_ECI[3][3]); 

void Dipole_mapping()

//Matrix operations 
void cross_prod(float *a, float *b, float *c);
void dot_prod(float *a, float *b, float *c);
void matmul(int m, int n, int p, float A[m][n], float B[n][p], float C[m][p]);
void vector_add(float *a, float *b, float *c);
void vector_times_scalar(float *a, float b, float *c);


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

//For calculationg the Quaternions 
float epsilon[3]={1,2,3}; float eta=1;  //These are the initial values of the quaternions
float eta_dot; float epsilon_dot[3]; 
float C_b_ECI[3][3];



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
    //Insert step 3 here. 
    quaternions(omega, epsilon, eta, &eta_dot, epsilon_dot);
    printf("dot product is %f %f %f %f", epsilon_dot[0], epsilon_dot[1], epsilon_dot[2], eta_dot);
	return (0);
}



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


void quaternions(float omega[], float epsilon[], float eta, float* eta_dot, float *epsilon_dot) {
	// 	This function takes in the initial values of eta and epsilon and calculates eta_dot and epsilon_dot
	// eta_dot and epsilon_dot then get integrated using RK4 to find eta and epsilon, repeat. 
	
	//Calculating eta
	float c_dot; 
	dot_prod(omega,epsilon, &c_dot);
	*eta_dot=-0.5*c_dot; 
	
	//Calculating Epsilon
	float epsilon_cross_omega[3];
	float eta_times_omega[3];
    float sum[3]; 

	cross_prod(epsilon, omega, epsilon_cross_omega);
    vector_times_scalar(omega, eta, eta_times_omega);
    vector_add(epsilon_cross_omega, eta_times_omega,sum);
    vector_times_scalar(sum, 0.5, epsilon_dot);	
}


void Rotation_Matrix_Body(float epsilon[], float eta, float C_b_ECI[3][3]){
	float e1=epsilon[0]; float e2=epsilon[1]; float e3=epsilon[2];
	//First row  
	C_b_ECI[1][1]=1-2*(e2*e2+e3*e3);
	C_b_ECI[1][2]=2*(e1*eta+e3*eta);
	C_b_ECI[1][3]=2*(e1*e3-e2*eta);
	
}


void Dipole_mapping(){
}



//Matrix Math functions 
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

void vector_add(float *a, float *b, float *c)
{
    c[0] = a[0]+b[0];
    c[1] = a[1]+b[1];
    c[2] = a[2]+b[2];
	}
	
void vector_times_scalar(float *a, float b, float *c)
{
	//a is the vector, b is the scalar 
    c[0] = a[0]*b;
    c[1] = a[1]*b;
    c[2] = a[2]*b;
	}