
/*** Matrix Math ***/
void matmul_scalar(int m, int n, float A[][n], float mult);
void matdiv_scalar(int m, int n, float A[][n], float div);
void transpose(int m, int n, float A[][n], float A_T[][m]);
void matmul(int m, int n, int p, float A[m][n], float B[n][p], float C[m][p]);
/*** Vector Math  ***/
void vector_add(float *a, float *b, float *c);
void vector_sub(float *a, float *b, float *c);
void vector_mult_scalar(float *v, float mult);
void vector_div_scalar(float *v, float div);
void dot_prod(float *a, float *b, float *c);
void cross_prod(float *a, float *b, float *c);
/*** Utilities ***/
void print_vector(int n, float *v);
void print_matrix(int m, int n, float A[][n]);
void solve(int m, float A[m][m], float b[m], float x[m]);
void invert(int m, float A[m][m], float AInv[m][m]);
void solveInv(int m, float A[m][m], float B[m][m], float X[m][m]);
