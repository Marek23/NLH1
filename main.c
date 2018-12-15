#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"

#define printfFnc(...) { mexPrintf(__VA_ARGS__); mexEvalString("drawnow;");}

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
            printfFnc("%.1f ", w[i][j]);
        }
    printfFnc("\n");
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
                printfFnc("%.1f ", w[i][j][k]);
            }
        printfFnc("\n");
        }
    printfFnc("\n \n");
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
                printfFnc("%d %d %d \n", m_it, n_it, k);
                for(i=0; i<s_s; i++)
                {
                    for(j=0; j<s_s; j++)
                    {
                        printfFnc("%f ", w[i][j][m_it][n_it][k]);
                    }
                    printfFnc("\n");
                }
                printfFnc("\n");
            }
            printfFnc("\n");
        }
        printfFnc("\n");
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

void initU(
    int m,
    int n,
    int c,
    int t_r,
    double in[m][n][c],
    double new_in[m+(2*t_r)][n+(2*t_r)][c])
    {
    printfFnc("INICJALIZACJA U \n");
    int new_m = m+(2*t_r);
    int new_n = n+(2*t_r);

    int i,j,k,i0,j0;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            i0=i+t_r;
            j0=j+t_r;
            for(k=0;k<c;k++)
            {
                new_in[i0][j0][k]=in[i][j][k];
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

int allOne(
    int m,
    int n,
    double in[m][n])
    {
    int i ,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(in[i][j]<1)
            {
                return 0;
            }
        }
    }
    return 1;
    }

void updatePhi(
    int M,
    int N,
    int c,
    double u[M][N][c],
    double phi[M][N])
    {
        int i,j;
        for(i=0;i<M;i++)
        {
            for(j=0;j<N;j++)
            {
                if(u[i][j][0]!=0   &&
                   u[i][j][1]!=255 &&
                   u[i][j][2]!=0   &&
                   phi[i][j]  ==0)
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
    float w[s_s][s_s][m][n][c],
    int iw,
    int jw,
    int kw,
    double out[s_s][s_s])
    {
        int i,j;
        for(i=0;i<s_s;i++)
        {
            for(j=0;j<s_s;j++)
            {
                out[i][j]=w[i][j][iw][jw][kw];
            }
        }
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
    float w[s_s][s_s][m][n][c])
    {
    printfFnc("PIERWSZE OBLICZENIE WAGI \n");
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
                                if (sw!=1)
                                {
                                    w[ii][jj][i][j][k]=
                                        exp((-1.0)*sumMatrix(2*p_r+1,2*p_r+1,ret )/(h*h))*
                                        exp((-1.0)*sumMatrix(2*k_r+1,2*k_r+1,retk)/(h*h));
                                }
                                else
                                {
                                    w[ii][jj][i][j][k]=
                                        exp((-1.0)*sumMatrix(2*p_r+1,2*p_r+1,ret )/(h*h));
                                }
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
    float w[s_s][s_s][m][n][c])
    {
    printfFnc("AKTUALIZACJA WAGI \n");
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

                if (PHI[i0][j0]==0)
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
                                    if (sw!=1)
                                    {
                                        w[ii][jj][i][j][k]=
                                            exp((-1.0)*sumMatrix(2*p_r+1,2*p_r+1,ret )/(h*h))*
                                            exp((-1.0)*sumMatrix(2*k_r+1,2*k_r+1,retk)/(h*h));
                                    }
                                    else
                                    {
                                        w[ii][jj][i][j][k]=
                                            exp((-1.0)*sumMatrix(2*p_r+1,2*p_r+1,ret )/(h*h));
                                    }
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

void solveNLCTV(
    int m,
    int n,
    int c,
    double u0[m][n][c],
    int M,
    int N,
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
    float w[s_s][s_s][m][n][c],
    double lamda,
    double f0[m][n][c])
    {
    int i,j,k,i0,j0,step;
    double u [M][N][c];

    initU(m,n,c,t_r,u0,u);

    for(step=1;step<99999;step++)
    {
        padarray3d(m,n,c,t_r,u0,u);
        if(step==1)
        {
            updateWeight(m,n,c,u0,M,N,u,h,p_s,kernel,
                p_sw,kernelk,t_r,s_r,p_r,sw,phi,s_s,w);
        }
        if(step>9 && step%10==0)
        {
            printfFnc("Step: %d \n", step);
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

                    if(phi[i0][j0]<1)
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
        if (allOne(M,N,phi)==1)
        {
            printfFnc("Koniec step: %d", step);
            break;
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
    int h            = mxGetScalar(pmxIn[6]);
    int p_s          = mxGetScalar(pmxIn[7]);
    double *kernelu  = mxGetPr(pmxIn[8]);
    int p_sw         = mxGetScalar(pmxIn[9]);
    double *kernelku = mxGetPr(pmxIn[10]);
    int t_r          = mxGetScalar(pmxIn[11]);
    int s_r          = mxGetScalar(pmxIn[12]);
    int p_r          = mxGetScalar(pmxIn[13]);
    int sw           = mxGetScalar(pmxIn[14]);
    double *phiu     = mxGetPr(pmxIn[15]);
    double *PHIu     = mxGetPr(pmxIn[16]);
    int s_s          = mxGetScalar(pmxIn[17]);
    double lamda     = mxGetScalar(pmxIn[18]);
    double *f0u      = mxGetPr(pmxIn[19]);

    int i,j,it,in4d,in2d,k,l,c_it,si,sj;

    double u0[m][n][c];
    double f0[m][n][c];

    double kernel [p_s] [p_s];
    double kernelk[p_sw][p_sw];

    double phi[M][N], PHI[M][N];

    printfFnc("UTWORZENIE ZMIENNEJ WAGI: ");
    float w[s_s][s_s][m][n][c];
    for(si=0;si<s_s;si++)
    {
        for(sj=0;sj<s_s;sj++)
        {
            for(i=0;i<m;i++)
            {
                for(j=0;j<n;j++)
                {
                    for(k=0;k<c;k++)
                    {
                        w[si][sj][i][j][k]=0;
                    }
                }
            }
        }
    }
    printfFnc("OK \n");
    c_it=0;
    it  =0;
    while(c_it<c)
    {
        while(it<(c_it+1)*m*n)
        {
            in2d=0;
            while(in2d<m*n)
            {
                u0[in2d%m][in2d/m][c_it]=u0u[it];
                f0[in2d%m][in2d/m][c_it]=f0u[it++];
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
        PHI[i%M][i/M] = PHIu[i];
    }

    solveNLCTV(m,n,c,u0,M,N,h,p_s,kernel,p_sw,kernelk,t_r,
               s_r,p_r,sw,phi,PHI,s_s,w,lamda,f0);

    pmxOut[0] = mxCreateDoubleMatrix(1,m*n*c,mxREAL);
    double *ret1;
    ret1 = mxGetPr(pmxOut[0]);

    c_it=0;
    it  =0;
    while(c_it<c)
    {
        while(it<(c_it+1)*m*n)
        {
            in2d=0;
            while(in2d<m*n)
            {
                ret1[it++] = u0[in2d%m][in2d/m][c_it];
                in2d++;
            }
        }
        c_it++;
    }
}
