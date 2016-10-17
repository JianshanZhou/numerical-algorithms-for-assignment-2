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
const int L = 100;//the tolerant iteration number


//declarations of all the sub-functions
int exchange_two_row(double ** array, int rn, int cn, int k, int ik);
int max_index(double **array, int rn, int cn, int k);
int Solve_linear_equations(double **A, double lambda, double *solution);
int Gauss_elimination(double **array, double *solution, int rn, int cn);

int main()
{
    printf("Hello world!\n");
    return 0;
}




/*This function initializes the matrix.*/
int initialize_matrix(double **A)
{
    return 0;
}


/*This function can be used to solve a specified linear equation system by
using the Gauss elimination algorithm. Note that the input array is the
corresponding augmented matrix.*/
int Solve_linear_equations(double **A, double lambda, double *solution)
{
    double array[N][N+1]={0.0};
    int rn=N, cn=(N+1);
    int i,j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<(N+1);j++)
        {
            if(j<N)
            {
                if(i==j)
                {
                    *((double*)array+cn*i+j)=*((double*)A+N*i+j)-lambda;
                }
                else
                {
                    *((double*)array+cn*i+j)=*((double*)A+N*i+j);
                }
            }
        }
    }
    Gauss_elimination(array, solution, rn, cn);
    return 0;
}


int Gauss_elimination(double **array, double *solution, int rn, int cn)
{

    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }

    int n=rn;

    /*step-1: eliminate elements of the given augmented matrix array*/
    int k, i, j, ak, ai, aj, ik;
    double m_ik=0.0;
    for(k=1; k<=(n-1); k++)
    {
        /*calculate the actual index k in C*/
        ak = k-1;
        /*find the index of the maximum element in the k-th column*/
        ik = max_index(array, rn, cn, k);
        /*exchange the k-th and the ik-th rows*/
        exchange_two_row(array, rn, cn, k, ik);

        for(i=(k+1); i<=n; i++)
        {
            /*calculate the actual index i in C*/
            ai = i-1;

            /*update the element values*/
            m_ik=(1.0*(*((double*)array+cn*ai+ak)))/(*((double*)array+cn*(ak)+ak));
            for(j=k+1; j<=n; j++)
            {
                /*calculate the actual index j in C*/
                aj = j-1;

                *((double*)array+cn*ai+aj) = *((double*)array+cn*ai+aj) - m_ik*(*((double*)array+cn*ak+aj));
            }
            *((double*)array+cn*ai+(n)) = *((double*)array+cn*ai+(n)) - m_ik*(*((double*)array+cn*ak+(n)));
        }
    }

    /*step-2: re-iteration*/
    double sum_temp=0.0;
    for(k=n; k>=1; k--)
    {
        ak = k-1;
        if(k==n)
        {
            *(solution + ak) = 1.0*(*((double*)array+cn*ak+n))/(*((double*)array+cn*ak+ak));
        }
        else
        {
            sum_temp = 0.0;
            for(j=k+1; j<=n; j++)
            {
                aj = j-1;
                sum_temp+=(*((double*)array+cn*ak+aj))*(*(solution + aj));
            }
            *(solution + ak) = 1.0*(*((double*)array+cn*ak+n)-sum_temp)/(*((double*)array+cn*ak+ak));
        }
    }

    printf("Successfully achieve the algorithm! The solution is:\n");
    for(i=1; i<=n; i++)
    {
        ai = i-1;
        printf("x[%d]=%lf\n",i,*(solution+ai));
    }
    return 0;
}


/*this function is used to find the index of the maximum element in a column vector:
k is the given index of the column of the given array*/
int max_index(double **array, int rn, int cn, int k)
{
    int i, ai, ak, ik;
    i = k;
    ai = i-1;
    ak = k-1;
    ik = i;

    double temp_max_element = fabs(*((double*)array+cn*ai+ak)), temp_value;
    for(i=k;i<=rn;i++){
        ai = i-1;
        temp_value = fabs(*((double*)array+cn*ai+ak));
        if(temp_max_element<temp_value){
            temp_max_element = temp_value;
            ik = i;
        }
    }
    return ik;/*output the index of the maximum element in the k-th column*/
}


/*this function is used to respectively exchange the corresponding elements in
the given k-th and the ik-th rows of the array.*/
int exchange_two_row(double ** array, int rn, int cn, int k, int ik)
{
    int ak = k-1, a_ik = ik-1, j, aj;
    int n=cn-1;
    double temp_value = 0.0;
    for(j=k;j<=rn;j++){
        aj = j-1;
        temp_value = *((double*)array + cn*ak + aj);
        *((double*)array + cn*ak + aj) =*((double*)array + cn*a_ik + aj);
        *((double*)array + cn*a_ik + aj) = temp_value;
    }
    temp_value = *((double*)array + cn*ak + n);
    *((double*)array + cn*ak + n) = *((double*)array + cn*a_ik + n);
    *((double*)array + cn*a_ik + n) = temp_value;
    return 0;
}
