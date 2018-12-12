#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include "mex.h"


void showMatrix2D(
    int m,
    int n,
    double w[m][n])
    {
    int i;
    int j;

    for(i=0; i<m; i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%.1f ", w[i][j]);
        }
    printf("\n");
    }
    }

void showMatrix3D(
    int m,
    int n,
    int c,
    double w[m][n][c])
    {
    int i;
    int j;
    int k;
    for(k=0; k<c; k++)
    {
        for(i=0; i<m; i++)
        {
            for(j=0;j<n;j++)
            {
                printf("%.1f ", w[i][j][k]);
            }
        printf("\n");
        }
    printf("\n \n");
    }
    }

void showMatrixWeight(
    int s_s,
    int m,
    int n,
    int c,
    double w[s_s][s_s][m][n][c])
    {
    int i;
    int j;
    int k;
    int m_it;
    int n_it;

    for(k=0; k<c; k++)
    {
        for(n_it=0; n_it<n; n_it++)
        {
            for(m_it=0; m_it<m; m_it++)
            {
                printf("%d %d %d \n", m_it, n_it, k);
                for(i=0; i<s_s; i++)
                {
                    for(j=0; j<s_s; j++)
                    {
                        printf("%f ", w[i][j][m_it][n_it][k]);
                    }
                    printf("\n");
                }
                printf("\n");
            }
            printf("\n");
        }
        printf("\n");
    }
    }
void padarray2d(
    int m,
    int n,
    int t_r,
    double in[m][n],
    double new_in[m+(2*t_r)][n+(2*t_r)])
    {
    int new_m = m+(2*t_r);
    int new_n = n+(2*t_r);

    int i;
    int j;

    int old_i;
    int old_j;

    for(old_i=0;old_i<m;old_i++)
    {
        for(old_j=0;old_j<n;old_j++)
        {
            i=old_i+t_r;
            j=old_j+t_r;

            new_in[i][j]=in[old_i][old_j];
        }
    }

    for(i=0;i<t_r;i++)
    {
        for(j=0;j<new_n;j++)
        {
            new_in[i][j]=new_in[2*t_r-i][j];
            new_in[new_m-t_r+i][j]=new_in[new_m-t_r-2-i][j];
        }
    }

    for(i=0;i<new_m;i++)
    {
        for(j=0;j<t_r;j++)
        {
            new_in[i][j]=new_in[i][2*t_r-j];
            new_in[i][new_n-t_r+j]=new_in[i][new_n-t_r-2-j];
        }
    }
    }

void padarray3d(
    int m,
    int n,
    int c,
    int t_r,
    double in[m][n][c],
    double new_in[m+(2*t_r)][n+(2*t_r)][c])
    {
    int new_m = m+(2*t_r);
    int new_n = n+(2*t_r);

    int i;
    int j;
    int k;
    int old_i;
    int old_j;

    for(old_i=0;old_i<m;old_i++)
    {
        for(old_j=0;old_j<n;old_j++)
        {
            i=old_i+t_r;
            j=old_j+t_r;
            for(k=0;k<c;k++)
            {
                new_in[i][j][k]=in[old_i][old_j][k];
            }
        }
    }

    for(i=0;i<t_r;i++)
    {
        for(j=0;j<new_n;j++)
        {
            for(k=0;k<c;k++)
            {
                new_in[i][j][k]=new_in[2*t_r-i][j][k];
                new_in[new_m-t_r+i][j][k]=new_in[new_m-t_r-2-i][j][k];
            }
        }
    }

    for(i=0;i<new_m;i++)
    {
        for(j=0;j<t_r;j++)
        {
            for(k=0;k<c;k++)
            {
                new_in[i][j][k]=new_in[i][2*t_r-j][k];
                new_in[i][new_n-t_r+j][k]=new_in[i][new_n-t_r-2-j][k];
            }
        }
    }
    }

void subMatrix3D(
    int p_r,
    int m,
    int n,
    int c,
    int k,
    double in[m][n][c],
    int i0,
    int j0,
    double out[2*p_r+1][2*p_r+1])
    {
    int i ,j;
    for(i=0;i<(2*p_r+1);i++)
    {
        for(j=0;j<(2*p_r+1);j++)
        {
            out[i][j]=in[i0-p_r+i][j0-p_r+j][k];
        }
    }
    }

void subMatrix(
    int p_r,
    int m,
    int n,
    double in[m][n],
    int i0,
    int j0,
    double out[2*p_r+1][2*p_r+1])
    {
    int i ,j;
    for(i=0;i<(2*p_r+1);i++)
    {
        for(j=0;j<(2*p_r+1);j++)
        {
            out[i][j]=in[i0-p_r+i][j0-p_r+j];
        }
    }
    }

void prodMatrix(
    int m,
    int n,
    double in1[m][n],
    double in2[m][n],
    double prod[m][n])
    {
    int i ,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            prod[i][j]=in1[i][j]*in2[i][j];
        }
    }
    }

void prod3Matrix(
    int m,
    int n,
    double in1[m][n],
    double in2[m][n],
    double in3[m][n],
    double prod[m][n])
    {
    int i ,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            prod[i][j]=in1[i][j]*in2[i][j]*in3[i][j];
        }
    }
    }

void subt2Matrix(
    int m,
    int n,
    double from[m][n],
    double q[m][n],
    double out[m][n])
    {
    int i ,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            out[i][j]=(from[i][j]-q[i][j])*(from[i][j]-q[i][j]);
        }
    }
    }

double sumMatrix(
    int m,
    int n,
    double in[m][n])
    {
    int i ,j;
    double sum=0;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            sum+=in[i][j];
        }
    }
    return sum;
    }

int any(
    int m,
    int n,
    double in[m][n])
    {
    int i ,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(in[i][j]>0)
            {
                return 1;
            }
        }
    }
    return 0;
    }

void updateWeight(
    int m,
    int n,
    int c,
    double u0[m][n][c],
    int M,
    int N,
    double u[M][N][c],
    int h,
    int p_s,
    double kernel[p_s][p_s],
    int p_sw,
    double kernelk[p_sw][p_sw],
    int t_r,
    int s_r,
    int p_r,
    int sw,
    double phi[M][N],
    int s_s,
    double w[s_s][s_s][m][n][c])
    {
    int k_r=sw*p_r;

    int i, j, k, r, s, i0, j0, ii, jj;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<c;k++)
            {
                i0=i+t_r;
                j0=j+t_r;
                double sphi [2*p_r+1][2*p_r+1];
                double sphik[2*k_r+1][2*k_r+1];
                subMatrix   (p_r,M,N,phi,i0,j0,sphi);
                subMatrix   (k_r,M,N,phi,i0,j0,sphik);

                if(any(2*p_r+1,2*p_r+1,sphi) == 1)
                {
                    double W1[2*p_r+1][2*p_r+1];
                    double W1k[2*k_r+1][2*k_r+1];
                    subMatrix3D(p_r,M,N,c,k,u,i0,j0,W1);
                    subMatrix3D(k_r,M,N,c,k,u,i0,j0,W1k);

                    ii=0;
                    for(r=i0-s_r;r<=i0+s_r;r++)
                    {
                        jj=0;
                        for(s=j0-s_r;s<=j0+s_r;s++)
                        {
                            if(phi[r][s]!=0)
                            {
                                double W2 [2*p_r+1][2*p_r+1];
                                double W2k[2*k_r+1][2*k_r+1];
                                subMatrix3D(p_r,M,N,c,k,u,r,s,W2);
                                subMatrix3D(k_r,M,N,c,k,u,r,s,W2k);

                                double diff [2*p_r+1][2*p_r+1];
                                double diffk[2*k_r+1][2*k_r+1];
                                subt2Matrix (2*p_r+1,2*p_r+1,W1,  W2  ,diff);
                                subt2Matrix (2*k_r+1,2*k_r+1,W1k, W2k ,diffk);

                                double ret [2*p_r+1][2*p_r+1];
                                double retk[2*k_r+1][2*k_r+1];
                                prod3Matrix(2*p_r+1,2*p_r+1,sphi ,kernel ,diff ,ret);
                                prod3Matrix(2*k_r+1,2*k_r+1,sphik,kernelk,diffk,retk);
                                w[ii][jj][i][j][k]=
                                    exp((-1.0)*sumMatrix(2*p_r+1,2*p_r+1,ret )/(h*h))*
                                    exp((-1.0)*sumMatrix(2*k_r+1,2*k_r+1,retk)/(h*h));
                            }
                            jj++;
                        }
                        ii++;
                    }
                }
            }
        }
    }
    }

void updateWeight2(
    int m,
    int n,
    int c,
    double u0[m][n][c],
    int M,
    int N,
    double u[M][N][c],
    int h,
    int p_s,
    double kernel[p_s][p_s],
    int p_sw,
    double kernelk[p_sw][p_sw],
    int t_r,
    int s_r,
    int p_r,
    int sw,
    double phi[M][N],
    double PHI[M][N],
    int s_s,
    double w[s_s][s_s][m][n][c])
    {
    int k_r=sw*p_r;

    int i, j, k, r, s, i0, j0, ii, jj;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<c;k++)
            {
                i0=i+t_r;
                j0=j+t_r;
                double sphi [2*p_r+1][2*p_r+1];
                double sphik[2*k_r+1][2*k_r+1];
                subMatrix   (p_r,M,N,phi,i0,j0,sphi);
                subMatrix   (k_r,M,N,phi,i0,j0,sphik);

                if (PHI(i0,j0)==0)
                {
                    if(any(2*p_r+1,2*p_r+1,sphi) == 1)
                    {
                        double W1[2*p_r+1][2*p_r+1];
                        double W1k[2*k_r+1][2*k_r+1];
                        subMatrix3D(p_r,M,N,c,k,u,i0,j0,W1);
                        subMatrix3D(k_r,M,N,c,k,u,i0,j0,W1k);

                        ii=0;
                        for(r=i0-s_r;r<=i0+s_r;r++)
                        {
                            jj=0;
                            for(s=j0-s_r;s<=j0+s_r;s++)
                            {
                                if(phi[r][s]!=0)
                                {
                                    double W2 [2*p_r+1][2*p_r+1];
                                    double W2k[2*k_r+1][2*k_r+1];
                                    subMatrix3D(p_r,M,N,c,k,u,r,s,W2);
                                    subMatrix3D(k_r,M,N,c,k,u,r,s,W2k);

                                    double diff [2*p_r+1][2*p_r+1];
                                    double diffk[2*k_r+1][2*k_r+1];
                                    subt2Matrix (2*p_r+1,2*p_r+1,W1,  W2  ,diff);
                                    subt2Matrix (2*k_r+1,2*k_r+1,W1k, W2k ,diffk);

                                    double ret [2*p_r+1][2*p_r+1];
                                    double retk[2*k_r+1][2*k_r+1];
                                    prod3Matrix(2*p_r+1,2*p_r+1,sphi ,kernel ,diff ,ret);
                                    prod3Matrix(2*k_r+1,2*k_r+1,sphik,kernelk,diffk,retk);
                                    w[ii][jj][i][j][k]=
                                        exp((-1.0)*sumMatrix(2*p_r+1,2*p_r+1,ret )/(h*h))*
                                        exp((-1.0)*sumMatrix(2*k_r+1,2*k_r+1,retk)/(h*h));
                                }
                                jj++;
                            }
                            ii++;
                        }
                    }
                }
            }
        }
    }
    }

void updatePhi(
    int M,
    int N,
    int c,
    double u[M][N][c],
    double phi[M][N])
    {
        int i,j,k,pom;
        for(i=0;i<M;i++)
        {
            for(j=0;j<N;j++)
            {
                if(u[i][j][0]!=0   ||
                   u[i][j][1]!=255 ||
                   u[i][j][2]!=0)
                {
                    phi[i][j]=1;
                }
            }
        }
    }

void subWeight(
    int m,
    int n,
    int c,
    int s_s,
    double w[s_s][s_s][m][n][c],
    int iw,
    int jw,
    int kw,
    double out[s_s][s_s])
    {
        int i,j
    for(i=0;i<s_s;i++)
    {
        for(j=0;j<s_s;j++)
        {
            out[i][j]=w[i][j][iw][jw][kw];
        }
    }
    }

void solveNLCTV(
    int step,
    int m,
    int n,
    int c,
    double u0[m][n][c],
    int M,
    int N,
    double u[M][N][c],
    int h,
    int p_s,
    double kernel[p_s][p_s],
    int p_sw,
    double kernelk[p_sw][p_sw],
    int t_r,
    int s_r,
    int p_r,
    int sw,
    double phi[M][N],
    double PHI[M][N],
    int s_s,
    double w[s_s][s_s][m][n][c],
    double lambda,
    double f0[m][n][c])
    {
    int i,j,k,i0,j0;
    padarray3d(m,n,c,t_r,u0,u);
    if(step == 1)
    {
        updateWeight(m,n,c,u0,M,N,u,h,p_s,kernel,
            p_sw,kernelk,t_r,s_r,p_r,sw,phi,s_s,w);
    }
    if(step%10==0)
    {
        updatePhi(M,N,c,u,phi);
        updateWeight2(m,n,c,u0,M,N,u,h,p_s,kernel,
            p_sw,kernelk,t_r,s_r,p_r,sw,phi,PHI,s_s,w);
    }
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<c;k++)
            {
                i0=i+t_r;
                j0=j+t_r;

                if(phi(i0,j0)<1)
                {
                    double sphi [2*p_r+1][2*p_r+1];
                    subMatrix   (p_r,M,N,phi,i0,j0,sphi);
                    if(any(2*p_r+1,2*p_r+1,sphi)==1)
                    {
                        double subw[s_s][s_s];
                        double subu[s_s][s_s];
                        subWeight(m,n,c,s_s,w,i,j,k,subw);
                        subMatrix3D(s_r,M,N,c,k,u,i0,j0,subu);
                        double prod[s_s][s_s];
                        prodMatrix(s_s,s_s,subw,subu,prod);
                        double sum_yw, sum_w;
                        sum_yw=sumMatrix(s_s,s_s,prod);
                        sum_w =sumMatrix(s_s,s_s,subw);

                        u0[i][j][k]=(lamda*PHI[i0][j0]*f0[i][j][k]+sum_yw)/(lamda*PHI[i0][j0]+sum_w);
                    }
                    else
                    {
                        u0[i][j][k]=u0[i][j][k];
                    }
                }
            }
        }
    }
    }

void mexFunction(int numOut, mxArray *pmxOut[],
                 int numIn, const mxArray *pmxIn[])
{
    int m            = mxGetScalar(pmxIn[0]);
    int n            = mxGetScalar(pmxIn[1]);
    int c            = mxGetScalar(pmxIn[2]);
    double *u0u      = mxGetPr(pmxIn[3]);
    int M            = mxGetScalar(pmxIn[4]);
    int N            = mxGetScalar(pmxIn[5]);
    double *uu       = mxGetPr(pmxIn[6]);
    int h            = mxGetScalar(pmxIn[7]);
    int p_s          = mxGetScalar(pmxIn[8]);
    double *kernelu  = mxGetPr(pmxIn[9]);
    int p_sw         = mxGetScalar(pmxIn[10]);
    double *kernelku = mxGetPr(pmxIn[11]);
    int t_r          = mxGetScalar(pmxIn[12]);
    int s_r          = mxGetScalar(pmxIn[13]);
    int p_r          = mxGetScalar(pmxIn[14]);
    int sw           = mxGetScalar(pmxIn[15]);
    double *phiu     = mxGetPr(pmxIn[16]);
    int s_s          = mxGetScalar(pmxIn[17]);
    double *wu       = mxGetPr(pmxIn[18]);

    int i, it, in4d, in2d, k, l, c_it;

    double u0[m][n][c];
    double u [M][N][c];

    double kernel [p_s] [p_s];
    double kernelk[p_sw][p_sw];

    double phi[M][N];

    double w[s_s][s_s][m][n][c];

    c_it=0;
    it  =0;
    while(c_it<c)
    {
        while(it<(c_it+1)*m*n)
        {
            in2d=0;
            while(in2d<m*n)
            {
                u0[in2d%m][in2d/m][c_it]=u0u[it++];
                in2d++;
            }
        }
        c_it++;
    }
    c_it=0;
    it  =0;
    while(c_it<c)
    {
        while(it<(c_it+1)*M*N)
        {
            in2d=0;
            while(in2d<M*N)
            {
                u[in2d%M][in2d/M][c_it]=uu[it++];
                in2d++;
            }
        }
        c_it++;
    }
    for(i=0;i<p_s*p_s;i++)
    {
        kernel[i%p_s][i/p_s] = kernelu[i];
    }
    for(i=0;i<p_sw*p_sw;i++)
    {
        kernelk[i%p_s][i/p_s] = kernelku[i];
    }
    for(i=0;i<M*N;i++)
    {
        phi[i%M][i/M] = phiu[i];
    }
    it   = 0;
    in4d = 0;
    c_it = 0;
    while (c_it<c)
    {
        while (it < (c_it+1)*s_s*s_s*m*n)
        {
            k = in4d%m;
            l = in4d/m;

            in2d=0;
            while (in2d < s_s*s_s)
            {
                w[in2d%s_s][in2d/s_s][k][l][c_it]=wu[it++];
                in2d++;
            }
            in4d++;
        }
        c_it++;
        in4d=0;
    }


    updateWeight(m,n,c,u0,M,N,u,h,p_s,kernel,
        p_sw,kernelk,t_r,s_r,p_r,sw,phi,s_s,w);

    pmxOut[0] = mxCreateDoubleMatrix(1,s_s*s_s*m*n*c,mxREAL);//num
    double *ret;
    ret = mxGetPr(pmxOut[0]);

    it   = 0;
    in4d = 0;
    c_it = 0;
    while (c_it<c)
    {
        while (it < (c_it+1)*s_s*s_s*m*n)
        {
            k = in4d%m;
            l = in4d/m;

            in2d=0;
            while (in2d < s_s*s_s)
            {
                ret[it++]=w[in2d%s_s][in2d/s_s][k][l][c_it];
                in2d++;
            }
            in4d++;
        }
        c_it++;
        in4d=0;
    }

}
/*
void mexFunction(int numOut, mxArray *pmxOut[],
                 int numIn, const mxArray *pmxIn[])
{
    int m            = mxGetScalar(pmxIn[0]);
    int n            = mxGetScalar(pmxIn[1]);
    int c            = mxGetScalar(pmxIn[2]);
    int pad          = mxGetScalar(pmxIn[3]);
    double *test1p   = mxGetPr(pmxIn[4]);
    double *test2p   = mxGetPr(pmxIn[5]);
    double *test3p   = mxGetPr(pmxIn[6]);
    int p_r          = mxGetScalar(pmxIn[7]);
    int i_0          = mxGetScalar(pmxIn[8]);
    int j_0          = mxGetScalar(pmxIn[9]);

    i_0--;
    j_0--;

    double toAny[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 1 } };

    double test1[m][n];
    double test2[m][n];
    double test3[m][n][c];

    int i, it, c_it, in2d;
    for(i=0;i<m*n;i++)
    {
        test1[i%m][i/m] = test1p[i];
        test2[i%m][i/m] = test2p[i];
    }

    c_it=0;
    it  =0;
    while(c_it<c)
    {
        while(it<(c_it+1)*m*n)
        {
            in2d=0;
            while(in2d<m*n)
            {
                test3[in2d%m][in2d/m][c_it]=test3p[it++];
                in2d++;
            }
        }
        c_it++;
    }

    int M = 2*pad+m;
    int N = 2*pad+n;
    double test1r[M][N];
    double test2r[M][N];
    double test3r[M][N][c];
    double test4r[2*p_r+1][2*p_r+1];
    double test5r[2*p_r+1][2*p_r+1];
    double test6r[2*p_r+1][2*p_r+1];
    double test7r[2*p_r+1][2*p_r+1];
    double test8r[2*p_r+1][2*p_r+1];

    pmxOut[0] = mxCreateDoubleMatrix(1,M*N,mxREAL);
    pmxOut[1] = mxCreateDoubleMatrix(1,M*N,mxREAL);
    pmxOut[2] = mxCreateDoubleMatrix(1,M*N*c,mxREAL);
    pmxOut[3] = mxCreateDoubleMatrix(1,(2*p_r+1)*(2*p_r+1),mxREAL);
    pmxOut[4] = mxCreateDoubleMatrix(1,(2*p_r+1)*(2*p_r+1),mxREAL);
    pmxOut[5] = mxCreateDoubleMatrix(1,(2*p_r+1)*(2*p_r+1),mxREAL);
    pmxOut[6] = mxCreateDoubleMatrix(1,(2*p_r+1)*(2*p_r+1),mxREAL);
    pmxOut[7] = mxCreateDoubleMatrix(1,(2*p_r+1)*(2*p_r+1),mxREAL);
    pmxOut[8] = mxCreateDoubleMatrix(1,1,mxREAL);
    pmxOut[9] = mxCreateNumericMatrix(1, 1,mxINT8_CLASS, 0);

    double *ret1 = mxGetPr(pmxOut[0]);
    double *ret2 = mxGetPr(pmxOut[1]);
    double *ret3 = mxGetPr(pmxOut[2]);
    double *ret4 = mxGetPr(pmxOut[3]);
    double *ret5 = mxGetPr(pmxOut[4]);
    double *ret6 = mxGetPr(pmxOut[5]);
    double *ret7 = mxGetPr(pmxOut[6]);
    double *ret8 = mxGetPr(pmxOut[7]);
    double *ret9 = mxGetPr(pmxOut[8]);
    int    *ret10= mxGetPr(pmxOut[9]);

    padarray2d(m,n,pad,test1,test1r);
    padarray2d(m,n,pad,test2,test2r);
    padarray3d(m,n,c,pad,test3,test3r);
    subMatrix3D(p_r,m,n,c,1,test3,i_0,j_0,test4r);
    subMatrix(p_r,m,n,test1,i_0,j_0,test5r);
    prodMatrix((2*p_r+1),(2*p_r+1),test5r,test5r,test6r);
    prod3Matrix((2*p_r+1),(2*p_r+1),test5r,test5r,test5r,test7r);
    subt2Matrix((2*p_r+1),(2*p_r+1),test4r,test5r,test8r);
    *ret9 = sumMatrix((2*p_r+1),(2*p_r+1),test8r);
    *ret10 = any(3,3,toAny);

    int k, l;
    it =0;
    while (it < M*N)
    {
        k = it%M;
        l = it/M;
        ret1[it] =test1r[k][l];
        ret2[it] =test2r[k][l];
        it++;
    }

    c_it=0;

    it  =0;
    in2d=0;
    while(c_it<c)
    {
        while(it<(c_it+1)*M*N)
        {
            in2d=0;
            while(in2d<M*N)
            {
                ret3[it++]=test3r[in2d%M][in2d/M][c_it];
                in2d++;
            }
        }
        c_it++;
    }

    it = 0;
    while (it < (2*p_r+1)*(2*p_r+1))
    {
        k = it%(2*p_r+1);
        l = it/(2*p_r+1);
        ret4[it]   =test4r[k][l];
        ret5[it]   =test5r[k][l];
        ret6[it]   =test6r[k][l];
        ret7[it]   =test7r[k][l];
        ret8[it++] =test8r[k][l];
    }
}
*/
