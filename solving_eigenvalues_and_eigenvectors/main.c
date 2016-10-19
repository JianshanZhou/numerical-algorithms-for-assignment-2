/* information
Copyright (C) Mon Oct 17 22:51:08 2016  Jianshan Zhou

Contact: zhoujianshan@buaa.edu.cn	jianshanzhou@foxmail.com

Website: <https://github.com/JianshanZhou>

This program is free software: you can redistribute
 it and/or modify it under the terms of
 the GNU General Public License as published
 by the Free Software Foundation,
 either version 3 of the License,
 or (at your option) any later version.

This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program.
 If not, see <http://www.gnu.org/licenses/>.

This is the main file where the main function and other sub-functions are provided
 to analyze the 2-st issue presented in the Numerical Analysis class. Specifically,
 the goal of this program is to solve for all the eigenvalues and eigenvectors of a specified
 matrix by using the QR decomposition algorithm and the Gauss elimination algorithm.
*/


//import some necessary basic libs
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


//define some macros
#define N 10


//tolerant accuracy
const double epsilon = 1.0e-12;//the accuracy
const int L = 1000;//the tolerant iteration number


//declarations of all the sub-functions
int exchange_two_row(double ** array, double *b, int rn, int cn, int k, int ik);
int max_index(double **array, int rn, int cn, int k);
int Solve_linear_equations(double **A, double lambda, double *solution);
int Gauss_elimination(double **array, double *solution, int rn, int cn);

int initialize_matrix(double **A);
double d_r(double **A_r, int r);
double c_r(double d_r, double a);
int sgn(double x);
double h_r(double c_r, double a);
int product_matrix_and_vector(double **matrix, double *vector, double *output, int n);
int product_Tmatrix_and_vector(double **matrix, double *vector, double *output, int n);
double inner_product(double *vector1, double *vector2, int n);
int outer_product(double *vector1, double *vector2, double **matrix, int n);
int subtraction_2matrix(double **matrix1, double **matrix2, double **ouput, int n);
int subtraction_3matrix(double **matrix1, double **matrix2, double **matrix3, double **ouput, int n);
int subtraction_vector(double *vector1, double *vector2, double c, double * output, int n);
int construct_u_r(double **A_r, double c_r, int r, int n, double * u_r);
int product_constant_and_vector(double * vector, double c, int n);
int update_A_r(double **A_r, double *u_r, double h_r, int r, int n);
double max_of_vector(double **A_r, int r, int n);
int approximated_triangulation(double **A, int n);
int test_on_triangulation(void);

int solve_for_two_roots(double **A, int m, double *s1R, double *s1I, double *s2R, double *s2I);
int product_2matrix(double **A, double **B, double **C, int m);
int construct_M(double **A, double **M, int m, double s, double t);
double d_r2(double **B, int r, int m);
int construct_u_r2(double **B, double *u_r, double c_r, int r, int m);
int update_B_r(double **B_r, double *u_r, double *v_r, int m);
int update_A_r2(double **A_r, double *u_r, double h_r, int m);
double max_of_vector2(double **B_r, int r, int m);
int perform_QR_on_A(double **A, int m);
int test_on_eigenvalues_solution(void);

int get_QR(double **A, double **Q);
int test_on_QR_decomposition(void);

int get_eigenvector(double **A, double lambda_j, int j, double **eigVector);
int get_eigenvectors(double ** orignal_A, double **A, double **eigVectors);
int test_on_eigenvectors_solution(void);

int inf(void)
{
    printf("*******************************************************************\n");
    printf("The whole program to solve for the eigenvalues of a matrix is developed by:\n");
    printf("Zhou Jianshan BY1613123\n");
    printf("Email: jianshanzhou@foxmail.com\n");
    printf("Web: https://github.com/JianshanZhou\n");
    printf("Please feel free to contact me if you have any questions!\n");
    printf("*******************************************************************\n\n");
    return 0;
}

int main()
{
    inf();

    //test_on_triangulation();
    //test_on_QR_decomposition();
//    test_on_eigenvalues_solution();
    test_on_eigenvectors_solution();

    return 0;
}


int test_on_eigenvectors_solution(void)
{
    double oA[N][N]={0.0};
    double A[N][N]={0.0};
    double eigVectors[N][N]={0.0};

    initialize_matrix(oA);
    initialize_matrix(A);

    approximated_triangulation(A,N);

    get_eigenvectors(oA, A, eigVectors);

    return 0;
}

int get_eigenvectors(double ** oA, double **A, double **eigVectors)
{
    double eigenR[N] = {0.0};
    double eigenI[N] = {0.0};
    solve_for_eigenvalues(A, eigenR, eigenI);

    int j, i;
    double lambda_j=0.0;
    for(j=0;j<N;j++)
    {
        if(fabs(*(eigenI+j))<=epsilon)
        {
            lambda_j = *(eigenR+j);
            get_eigenvector(oA, lambda_j, j, eigVectors);

            printf("***The eigenvector of lambda[%d]=(%.11e,\t%.11e)***\n",j+1,*(eigenR+j),*(eigenI+j));
            for(i=0;i<N;i++)
            {
                printf("%.11e\n",*((double*)eigVectors+N*i+j));
            }
            printf("\n");
        }
    }


    return 0;
}

int get_eigenvector(double **A, double lambda_j, int j, double **eigVector)
{
    double solution[N] = {0.0};
    Solve_linear_equations(A, lambda_j, solution);
    int i;
    for(i=0;i<N;i++)
    {
        *((double*)eigVector+N*i+j) = *(solution+i);
    }
    return 0;
}


int test_on_QR_decomposition(void)
{
    double R[N][N] = {0.0};
    double Q[N][N] = {0.0};
    initialize_matrix(R);
    approximated_triangulation(R,N);

    get_QR(R,Q);

    int i, j;
    printf("****************The Q matrix*****************\n");
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("Q[%d][%d]=%.11e\t",i+1,j+1,Q[i][j]);
            if(((j+1)%2)==0)
            {
                printf("\n");
            }
        }
        printf("\n\n");
    }

    printf("****************The R matrix*****************\n");
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("R[%d][%d]=%.11e\t",i+1,j+1,R[i][j]);
            if(((j+1)%2)==0)
            {
                printf("\n");
            }
        }
        printf("\n\n");
    }

    return 0;
}


int get_QR(double **A, double **Q)
{
    //initialize Q as I
    int i, j, r;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if(i == j)
            {
                *((double*)Q+N*i+j) = 1.0;
            }
        }
    }

    //loop
    double maxV = 0.0, dd_r, cc_r, hh_r;
    double u_r[N]={0.0}, omega[N]={0.0}, p_r[N]={0.0};

    for(r=1;r<(N-1);r++)
    {
        maxV = max_of_vector2(A, r, N);
        if(fabs(maxV)<=epsilon)
        {
            ;
        }
        else
        {
            dd_r = d_r2(A, r, N);
            cc_r = c_r(dd_r,*((double*)A+N*(r-1)+(r-1)));
            hh_r = h_r(cc_r,*((double*)A+N*(r-1)+(r-1)));
            construct_u_r2(A,u_r,cc_r,r,N);

            product_matrix_and_vector(Q,u_r,omega,N);
            for(i=0;i<N;i++)
            {
                for(j=0;j<N;j++)
                {
                    *((double*)Q+N*i+j)=*((double*)Q+N*i+j)-(*(omega+i))*(*(u_r+j))/hh_r;
                }
            }

            product_Tmatrix_and_vector(A,u_r,p_r,N);
            product_constant_and_vector(p_r,1.0/hh_r,N);

            for(i=0;i<N;i++)
            {
                for(j=0;j<N;j++)
                {
                    *((double*)A+N*i+j)=*((double*)A+N*i+j)-(*(u_r+i))*(*(p_r+j));
                }
            }

        }
    }

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if(fabs(*((double*)Q+N*i+j))<=epsilon)
            {
                (*((double*)Q+N*i+j))=0.0;
            }
            if(fabs(*((double*)A+N*i+j))<=epsilon)
            {
                (*((double*)A+N*i+j))=0.0;
            }
        }
    }
    return 0;
}


int test_on_eigenvalues_solution(void)
{

    double A[N][N] = {0.0};
    initialize_matrix(A);

    approximated_triangulation(A, N);

    double eigenR[N] = {0.0};
    double eigenI[N] = {0.0};

    solve_for_eigenvalues(A, eigenR, eigenI);

    printf("*************************The final QR matrix*****************************\n");
    int i, j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("A[%d][%d]=%.11e\t",i+1,j+1,A[i][j]);
            if(((j+1)%2)==0)
            {
                printf("\n");
            }
        }
        printf("\n\n");
    }
    printf("*************************The eigenvalues*****************************\n\n");
    for(i=0;i<N;i++)
    {
        printf("lambda[%d]=(%.11e,\t%.11ei)\n",i+1,eigenR[i],eigenI[i]);
    }

    return 0;
}

/*This function solves for the eigenvalues of A based on two-shift QR decomposition algorithm.*/
int solve_for_eigenvalues(double **A, double *eigR, double *eigI)
{
    double s1R, s1I, s2R, s2I;
    int k=1, m=N;
    while((k<=L)&&(m>=1))
    {
        if(m == 1)//(4)
        {
            *(eigR + m-1) = *((double*)A+N*(m-1)+(m-1));
            *(eigI + m-1) = 0.0;
            m = m-1;
        }
        else
        {
            //(3)
            if(fabs(*((double*)A+N*(m-1)+(m-2)))<=epsilon)
            {
                *(eigR + m-1) = *((double*)A+N*(m-1)+(m-1));
                *(eigI + m-1) = 0.0;
                m = m-1;
            }
            else
            {
                //(5)
                solve_for_two_roots(A,m,&s1R,&s1I,&s2R,&s2I);
                //(6)
                if(m == 2)
                {
                    //solve for two roots

                    *(eigR + m-1-1) = s1R;
                    *(eigI + m-1-1) = s1I;
                    *(eigR + m-1) = s2R;
                    *(eigI + m-1) = s2I;
                    m = m-2;
                }
                else if(fabs(*((double*)A+N*(m-2)+(m-3)))<=epsilon)//(7)
                {
                    *(eigR + m-1-1) = s1R;
                    *(eigI + m-1-1) = s1I;
                    *(eigR + m-1) = s2R;
                    *(eigI + m-1) = s2I;
                    m = m-2;
                }
                else
                {
                    perform_QR_on_A(A,m);
                }
            }
        }
        k++;
    }

    return 0;
}


int perform_QR_on_A(double **A, int m)
{
    double s = (*((double*)A+N*(m-2)+(m-2)))+(*((double*)A+N*(m-1)+(m-1)));
    double t1 = (*((double*)A+N*(m-2)+(m-2)))*(*((double*)A+N*(m-1)+(m-1)));
    double t2 = (*((double*)A+N*(m-1)+(m-2)))*(*((double*)A+N*(m-2)+(m-1)));
    double t = t1-t2;

    double Mk[N][N] = {0.0};
    double u_r[N]={0.0}, v_r[N]={0.0};
    double maxV=0.0, dd_r, cc_r, hh_r;

    construct_M(A,Mk,m,s,t);

    int r;
    for(r=1;r<=(m-1);r++)
    {
        maxV = max_of_vector2(Mk, r, m);
        if(maxV<=epsilon)
        {
            ;
        }
        else
        {
            dd_r = d_r2(Mk,r,m);
            cc_r = c_r(dd_r,(*((double*)Mk+N*(r-1)+(r-1))));
            hh_r = h_r(cc_r, (*((double*)Mk+N*(r-1)+(r-1))));
            construct_u_r2(Mk,u_r,cc_r,r,m);
            product_Tmatrix_and_vector(Mk,u_r,v_r,m);
            product_constant_and_vector(v_r, 1.0/hh_r, m);
            update_B_r(Mk,u_r,v_r,m);
            update_A_r2(A,u_r,hh_r,m);
        }
    }
    return 0;
}

double max_of_vector2(double **B_r, int r, int m)
{
    double maxV=0.0;
    double temp;
    int i;
    for(i=(r+1);i<=m;i++)
    {
        temp = fabs(*((double*)B_r+N*(i-1)+(r-1)));
        if(temp>maxV)
        {
            maxV = temp;
        }

    }
    return maxV;
}

int update_A_r2(double **A_r, double *u_r, double h_r, int m)
{
    double p_r[N]={0.0};
    double q_r[N]={0.0};
    double omega[N]={0.0};
    double t_r=0.0;

    product_Tmatrix_and_vector(A_r,u_r,p_r,m);
    product_constant_and_vector(p_r,1.0/h_r,m);

    product_matrix_and_vector(A_r,u_r,q_r,m);
    product_constant_and_vector(q_r,1.0/h_r,m);

    t_r = inner_product(p_r,u_r,m)/h_r;

    int i,j;

    for(i=0;i<m;i++)
    {
        *(omega+i) = *(q_r+i) - t_r*(*(u_r+i));
    }

    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            *((double*)A_r+N*i+j) = *((double*)A_r+N*i+j) - (*(omega+i))*(*(u_r+j)) - (*(u_r+i))*(*(p_r+j));
        }
    }
    return 0;
}

int update_B_r(double **B_r, double *u_r, double *v_r, int m)
{
    int i,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            *((double*)B_r+N*i+j) = *((double*)B_r+N*i+j) - (*(u_r+i))*(*(v_r+j));
        }
    }
    return 0;
}

int construct_u_r2(double **B, double *u_r, double c_r, int r, int m)
{
    int i;
    for(i=1;i<=m;i++)
    {
        if(i<r)
        {
            *(u_r+i-1) = 0.0;
        }
        else if(i==r)
        {
            *(u_r+i-1) = *((double*)B+N*(i-1)+(r-1))-c_r;
        }
        else
        {
            *(u_r+i-1) = *((double*)B+N*(i-1)+(r-1));
        }
    }
    return 0;
}

double d_r2(double **B, int r, int m)
{
    double d=0.0;
    int i;
    for(i=r;i<=m;i++)
    {
        d += (*((double*)B+N*(i-1)+(r-1)))*(*((double*)B+N*(i-1)+(r-1)));
    }
    d = sqrt(d);
    return d;
}

int construct_M(double **A, double **M, int m, double s, double t)
{
    int i, j;
    double AA[N][N] = {0.0};
    product_2matrix(A, A, AA, m);
    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            if(i==j)
            {
                *((double*)M+N*i+j) = *((double*)AA+N*i+j) - s*(*((double*)A+N*i+j)) + t*1.0;
            }
            else
            {
                *((double*)M+N*i+j) = *((double*)AA+N*i+j) - s*(*((double*)A+N*i+j));
            }
        }
    }

    return 0;
}


int product_2matrix(double **A, double **B, double **C, int m)
{
    int i, j, t;
    double temp = 0.0;

    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            temp = 0.0;
            for(t=0;t<m;t++)
            {
                temp += (*((double*)A+N*i+t))*(*((double*)B+N*t+j));
            }
            *((double*)C+N*i+j) = temp;
        }
    }
    return 0;
}

int solve_for_two_roots(double **A, int m, double *s1R, double *s1I, double *s2R, double *s2I)
{

    double a11 = *((double*)A+N*(m-2)+(m-2));
    double a12 = *((double*)A+N*(m-2)+(m-1));
    double a21 = *((double*)A+N*(m-1)+(m-2));
    double a22 = *((double*)A+N*(m-1)+(m-1));

    double a = 1.0, b = -1.0*(a11+a22), c = a11*a22 - a12*a21;

    double delta = b*b - 4.0*a*c;

    if(delta<0)
    {
        *s1R = -1.0*b/(2.0*a);
        *s2R = -1.0*b/(2.0*a);
        *s1I = sqrt(-1.0*delta)/(2.0*a);
        *s2I = -1.0*sqrt(-1.0*delta)/(2.0*a);
    }
    else
    {
        *s1R = (-1.0*b+sqrt(delta))/(2.0*a);
        *s1I = 0.0;
        *s2R = (-1.0*b-sqrt(delta))/(2.0*a);
        *s2I = 0.0;
    }
    return 0;
}


/*a test function for approximative triangulation of the given matrix.*/
int test_on_triangulation(void)
{
    double A[N][N] = {0.0};
    initialize_matrix(A);

    approximated_triangulation(A, N);

    int i, j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("A[%d][%d]=%.11e\t",i+1,j+1,A[i][j]);
            if(((j+1)%2) == 0)
            {
                printf("\n");
            }
        }
        printf("\n\n");
    }
    return 0;
}

/*The following functions are used to approximatively triangulate the given matrix A.*/
int approximated_triangulation(double **A, int n)
{
    int r;
    double max_value = 0.0;
    double rh_r = 0.0, rc_r = 0.0;
    double u_r[N] = {0.0};

    for(r=1;r<=(n-2);r++)
    {
        max_value = max_of_vector(A, r, n);
        if(fabs(max_value)<=epsilon)
        {
            ;
        }
        else
        {
            rc_r = c_r(d_r(A, r), (*((double*)A+n*(r)+r-1)));
            rh_r = h_r(rc_r, (*((double*)A+n*r+r-1)));
            construct_u_r(A,rc_r,r,n,u_r);
            update_A_r(A,u_r,rh_r,r,n);
        }
    }



    return 0;
}

double max_of_vector(double **A_r, int r, int n)
{
    double temp = 0.0;
    double max_value = 0.0;
    int i;
    for(i=(r+2);i<=n;i++)
    {
        temp = fabs(*((double*)A_r+N*(i-1)+r-1));
        if(temp>max_value)
        {
            max_value = temp;
        }
    }
    return max_value;
}

int update_A_r(double **A_r, double *u_r, double h_r, int r, int n)
{
    double p_r[N] = {0.0};
    double q_r[N] = {0.0};
    double t_r = 0.0;
    double omega_r[N] = {0.0};

    //calculate p_r = A_r^T*u_r/h_r
    product_Tmatrix_and_vector(A_r, u_r, p_r, n);
    product_constant_and_vector(p_r, 1.0/h_r, n);
    //calculate q_r = A_r*u_r/h_r
    product_matrix_and_vector(A_r, u_r, q_r, n);
    product_constant_and_vector(q_r, 1.0/h_r, n);
    //calculate t_r
    t_r = inner_product(p_r, u_r, n)/h_r;
    //calculate omega_r
    subtraction_vector(q_r, u_r, t_r, omega_r, n);
    //omega_ru_r^T
    double omega_u[N][N] = {0.0};
    double u_p[N][N] = {0.0};
    outer_product(omega_r, u_r, omega_u, n);
    outer_product(u_r, p_r, u_p, n);

    subtraction_3matrix(A_r, omega_u, u_p, A_r, n);

    return 0;
}

int product_constant_and_vector(double * vector, double c, int n)
{
    int i;
    for(i=0;i<n;i++)
    {
        *(vector+i) = (c)*(*(vector+i));
    }
    return 0;
}

int construct_u_r(double **A_r, double c_r, int r, int n, double *u_r)
{
    int aj = r-1;
    int i;
    for(i=0;i<n;i++)
    {
        if(i<(r))
        {
            *(u_r + i) = 0.0;
        }
        else if(i == r)
        {
            *(u_r + i) = *((double*)A_r+N*i+aj)-c_r;
        }
        else
        {
            *(u_r + i) = *((double*)A_r+N*i+aj);
        }
    }
    return 0;
}

int subtraction_vector(double *vector1, double *vector2, double c, double *output, int n)
{
    int i;
    for(i=0;i<n;i++)
    {
        *(output + i) = *(vector1 + i) - (c)*(*(vector2 + i));
    }
    return 0;
}

int subtraction_3matrix(double **matrix1, double **matrix2, double **matrix3, double **ouput, int n)
{
    subtraction_2matrix(matrix1, matrix2, ouput, n);
    subtraction_2matrix(ouput, matrix3, ouput, n);
    return 0;
}

int subtraction_2matrix(double **matrix1, double **matrix2, double **ouput, int n)
{
    int i, j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            *((double*)ouput+N*i+j) = *((double*)matrix1+N*i+j) - *((double*)matrix2+N*i+j);
        }
    }
    return 0;
}

int outer_product(double *vector1, double *vector2, double **matrix, int n)
{
    int i, j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            *((double*)matrix + i*N + j) = (*(vector1 + i))*(*(vector2 + j));
        }
    }
    return 0;
}

double inner_product(double *vector1, double *vector2, int n)
{
    double temp = 0.0;
    int i;
    for(i=0;i<n;i++)
    {
        temp += (*(vector1 + i))*(*(vector2+i));
    }
    return temp;
}

double d_r(double **A_r, int r)
{
    double temp = 0.0;
    int i, ai;
    int aj = r-1;
    for(i=(r+1);i<=N;i++)
    {
        ai = i-1;
        temp += (*((double*)A_r+N*ai+aj))*(*((double*)A_r+N*ai+aj));
    }
    temp = sqrt(temp);
    return temp;
}


double c_r(double d_r, double a)
{
    double c;
    if(fabs(a)<=epsilon)
    {
        return d_r;
    }
    else
    {
        c = (-1.0)*sgn(a)*d_r;
        return c;
    }
}

int sgn(double x)
{
    int s;
    if(x>0)
    {
        s = 1;
    }
    else if(x<0)
    {
        s = -1;
    }
    else
    {
        s = 0;
    }
    return s;
}

double h_r(double c_r, double a)
{
    return c_r*c_r - c_r*a;
}

int product_matrix_and_vector(double **matrix, double *vector, double *output, int n)
{
    int i, j;
    double temp = 0.0;
    for(i=0;i<n;i++)
    {
        temp = 0.0;
        for(j=0;j<n;j++)
        {
            temp += (*((double*)matrix + N*i +j))*(*(vector + j));
        }
        *(output+i) = temp;
    }

    return 0;
}

int product_Tmatrix_and_vector(double **matrix, double *vector, double *output, int n)
{
    int i, j;
    double temp = 0.0;
    for(j=0;j<n;j++)
    {
        temp = 0.0;
        for(i=0;i<n;i++)
        {
            temp += (*((double*)matrix + N*i + j))*(*(vector+i));
        }
        *(output+j) = temp;
    }
    return 0;
}


/*This function initializes the matrix.*/
int initialize_matrix(double **A)
{
    int i, j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if(i == j)
            {
                *((double*)A+N*i+j) = 1.52*cos(1.0*(i+1)+1.2*(j+1));
            }
            else
            {
                *((double*)A+N*i+j) = sin(0.5*(i+1)+0.2*(j+1));
            }
        }
    }
    return 0;
}



/*This function can be used to solve a specified linear equation system by
using the Gauss elimination algorithm. Note that the input array is the
corresponding augmented matrix.*/
int Solve_linear_equations(double **A, double lambda, double *solution)
{
    double array[N][N]={0.0};
    int i,j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            array[i][j] = *((double*)A+N*i+j);
            if(i==j)
            {
                array[i][j] = array[i][j] - lambda;
            }
        }
    }
    Gauss_elimination(array, solution, N, N);
    return 0;
}


int Gauss_elimination(double **array, double *solution, int rn, int cn)
{
    double b[N]={0.0};
    double mik;
    int k,ik,i,j;
    for(k=1;k<=(N-1);k++)
    {
        ik = max_index(array,rn,cn,k);
        exchange_two_row(array, b, rn, cn, k, ik);
        for(i=(k+1);i<=rn;i++)
        {
            mik = (*((double*)array+cn*(i-1)+(k-1)))/(*((double*)array+cn*(k-1)+(k-1)));
            for(j=(k+1);j<=cn;j++)
            {
                *((double*)array+cn*(i-1)+(j-1))=*((double*)array+cn*(i-1)+(j-1))-mik*(*((double*)array+cn*(k-1)+(j-1)));
            }
            *(b+i-1) = *(b+i-1)-mik*(*(b+k-1));
        }
    }

    if(fabs(*(b+rn-1))<=epsilon)
    {
        *(solution+rn-1) = 1.0;
    }
    else
    {
        *(solution+rn-1) = (*(b+rn-1))/(*((double*)array+cn*(rn-1)+(cn-1)));
    }

    double temp;
    for(k=(rn-1);k>=1;k--)
    {
        temp = 0.0;
        for(j=(k+1);j<=cn;j++)
        {
            temp += (*((double*)array+cn*(k-1)+(j-1)))*(*(solution+j-1));
        }
        *(solution+k-1) = (*(b+k-1)-temp)/(*((double*)array+cn*(k-1)+(k-1)));
    }

    return 0;
}


/*this function is used to find the index of the maximum element in a column vector:
k is the given index of the column of the given array*/
int max_index(double **array, int rn, int cn, int k)
{
    double maxV;
    int i=k,ik;
    maxV = fabs(*((double*)array+cn*(i-1)+(k-1)));
    ik = k;
    for(i=k;i<=rn;i++)
    {
        if(fabs(*((double*)array+cn*(i-1)+(k-1)))>maxV)
        {
            ik=i;
            maxV = fabs(*((double*)array+cn*(i-1)+(k-1)));
        }
    }
    return ik;/*output the index of the maximum element in the k-th column*/
}


/*this function is used to respectively exchange the corresponding elements in
the given k-th and the ik-th rows of the array.*/
int exchange_two_row(double ** array, double *b, int rn, int cn, int k, int ik)
{
    int j;
    double temp1, temp2;
    for(j=k;j<=cn;j++)
    {
        temp1 = *((double*)array+cn*(k-1)+(j-1));
        *((double*)array+cn*(k-1)+(j-1)) = *((double*)array+cn*(ik-1)+(j-1));
        *((double*)array+cn*(ik-1)+(j-1)) = temp1;
    }
    temp2 = *(b+k-1);
    *(b+k-1) = *(b+ik-1);
    *(b+ik-1) = temp2;
    return 0;
}
