#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"

#define printfFnc(...) { mexPrintf(__VA_ARGS__); mexEvalString("drawnow;");}

float **sphi,  **W1,  **W2,  **diff,  **ret,  **kernel;  // p_s x p_s
float **sphik, **W1k, **W2k, **diffk, **retk, **kernelk; // P_S x P_S

float **subw, **subu, **prod;                 // s_s x s_s

float **phi2, **phi, **PHI;                   // M x N

float ***u;                                   // M x N x c

float ***u0, ***f0;                           // m x n x c

float *****w;                                 // s_s x s_s x m x n

void initVars(int p_s, int P_S, int s_s, int M, int N, int m, int n, int c) {
    int iii, jjj, si, sj;
    sphi   = calloc(p_s, sizeof(float *));
    W1     = calloc(p_s, sizeof(float *));
    W2     = calloc(p_s, sizeof(float *));
    diff   = calloc(p_s, sizeof(float *));
    ret    = calloc(p_s, sizeof(float *));
    kernel = calloc(p_s, sizeof(float *));
    for(iii = 0; iii < p_s; iii++) {
        sphi[iii]   = calloc(p_s, sizeof(float));
        W1[iii]     = calloc(p_s, sizeof(float));
        W2[iii]     = calloc(p_s, sizeof(float));
        diff[iii]   = calloc(p_s, sizeof(float));
        ret[iii]    = calloc(p_s, sizeof(float));
        kernel[iii] = calloc(p_s, sizeof(float));
    }

    sphik   = calloc(P_S, sizeof(float *));
    W1k     = calloc(P_S, sizeof(float *));
    W2k     = calloc(P_S, sizeof(float *));
    diffk   = calloc(P_S, sizeof(float *));
    retk    = calloc(P_S, sizeof(float *));
    kernelk = calloc(P_S, sizeof(float *));
    for(iii = 0; iii < P_S; iii++) {
        sphik[iii]   = calloc(P_S, sizeof(float));
        W1k[iii]     = calloc(P_S, sizeof(float));
        W2k[iii]     = calloc(P_S, sizeof(float));
        diffk[iii]   = calloc(P_S, sizeof(float));
        retk[iii]    = calloc(P_S, sizeof(float));
        kernelk[iii] = calloc(P_S, sizeof(float));
    }

    subw = calloc(s_s, sizeof(float *));
    subu = calloc(s_s, sizeof(float *));
    prod = calloc(s_s, sizeof(float *));
    for(iii = 0; iii < s_s; iii++) {
        subw[iii] = calloc(s_s, sizeof(float));
        subu[iii] = calloc(s_s, sizeof(float));
        prod[iii] = calloc(s_s, sizeof(float));
    }

    phi2 = calloc(M, sizeof(float *));
    phi  = calloc(M, sizeof(float *));
    PHI  = calloc(M, sizeof(float *));
    for(iii = 0; iii < M; iii++) {
       phi2[iii] = calloc(N, sizeof(float));
       phi[iii]  = calloc(N, sizeof(float));
       PHI[iii]  = calloc(N, sizeof(float));
    }

    u = calloc(M, sizeof(float **));
    for(iii = 0; iii < M; iii++) { 
       u[iii] = calloc(N, sizeof(float *));
        for(jjj = 0; jjj < N; jjj++) { 
            u[iii][jjj] = calloc(c, sizeof(float));
        }
    }

    u0 = calloc(m, sizeof(float **));
    f0 = calloc(m, sizeof(float **));
    for(iii = 0; iii < m; iii++) { 
       u0[iii] = calloc(n, sizeof(float *));
       f0[iii] = calloc(n, sizeof(float *));
        for(jjj = 0; jjj < n; jjj++) { 
            u0[iii][jjj] = calloc(c, sizeof(float));
            f0[iii][jjj] = calloc(c, sizeof(float));
        }
    }

    printfFnc("Inicjalizacja wagi: ");
    w = calloc(s_s, sizeof(float ****));
    if(!w){
        printfFnc("FAILED \n");
        return;
    }
    for(si = 0; si < s_s; si++) { 
        w[si] = calloc(s_s, sizeof(float ***));
        if(!w[si]){
            printfFnc("FAILED \n");
            return;
        }
        for(sj = 0; sj < s_s; sj++) {
            w[si][sj] = calloc(m, sizeof(float **));
            if(!w[si][sj]){
                printfFnc("FAILED \n");
                return;
            }
            for(iii = 0; iii < m; iii++) { 
                w[si][sj][iii] = calloc(n, sizeof(float *));
                if(!w[si][sj][iii]){
                    printfFnc("FAILED \n");
                    return;
                }
                for(jjj = 0; jjj < n; jjj++) { 
                    w[si][sj][iii][jjj] = calloc(c, sizeof(float));
                    if(!w[si][sj][iii][jjj]){
                        printfFnc("FAILED \n");
                        return;
                    }
                }
            }
        }
    }
    printfFnc("OK. \n");
}

void clearVars(int p_s, int P_S, int s_s, int M, int N, int m, int n, int c) {
    int iii, jjj, si, sj;

    // float **sphi,  **W1,  **W2,  **diff,  **ret, **kernel;  // p_s x p_s
    for(iii = 0; iii < p_s; iii++) {
        free(sphi[iii]);
        free(W1[iii]);
        free(W2[iii]);
        free(diff[iii]);
        free(ret[iii]);
        free(kernel[iii]);
    }
    free(sphi);
    free(W1);
    free(W2);
    free(diff);
    free(ret);
    free(kernel);

    // float **sphik, **W1k, **W2k, **diffk, **retk, **kernelk; // P_S x P_S
    for(iii = 0; iii < P_S; iii++) {
        free(sphik[iii]);
        free(W1k[iii]);
        free(W2k[iii]);
        free(diffk[iii]);
        free(retk[iii]);
        free(kernelk[iii]);
    }
    free(sphik);
    free(W1k);
    free(W2k);
    free(diffk);
    free(retk);
    free(kernelk);

    // float **subw, **subu, **prod;                 // s_s x s_s
    for(iii = 0; iii < s_s; iii++) {
        free(subw[iii]);
        free(subu[iii]);
        free(prod[iii]);
    }
    free(subw);
    free(subu);
    free(prod);

    // float **phi2, **phi, **PHI;                   // M x N
    for(iii = 0; iii < M; iii++) {
        free(phi2[iii]);
        free(phi[iii]);
        free(PHI[iii]);
    }
    free(phi2);
    free(phi);
    free(PHI);

    // float ***u; // M x N x c
    for(iii = 0; iii < M; iii++) {
        for(jjj = 0; jjj < N; jjj++) {
            free(u[iii][jjj]);
        }
        free(u[iii]);
    }
    free(u);

    // float ***u0, ***f0;                           // m x n x c
    for(iii = 0; iii < m; iii++) {
        for(jjj = 0; jjj < n; jjj++) {
            free(u0[iii][jjj]);
            free(f0[iii][jjj]);
        }
        free(u0[iii]);
        free(f0[iii]);
    }
    free(u0);
    free(f0);

    // float *****w;                                 // s_s x s_s x m x n
    for(si = 0; si < s_s; si++) {
        for(sj = 0; sj < s_s; sj++) {
            for(iii = 0; iii < m; iii++) {
                for(jjj = 0; jjj < n; jjj++) {
                    free(w[si][sj][iii][jjj]);
                }
                free(w[si][sj][iii]);
            }
            free(w[si][sj]);
        }
        free(w[si]);
    }
    free(w);

}

void showMatrix2D(
    int m,
    int n,
    float **w)
    {
    int i;
    int j;

    for(i=0; i<m; i++)
    {
        for(j=0;j<n;j++)
        {
            printfFnc("%.1f ", w[i][j]);
        }
    printf("\n");
    }
    }

void showMatrix3D(
    int m,
    int n,
    int c,
    float ***w)
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
    float *****w)
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
    float **in,
    float **new_in)
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
    float ***in,
    float ***new_in)
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
    float ***in,
    int i0,
    int j0,
    float **out)
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
    float **in,
    int i0,
    int j0,
    float **out)
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
    float **in1,
    float **in2,
    float **prod)
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
    float **in1,
    float **in2,
    float **in3,
    float **prod)
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
    float **from,
    float **q,
    float **out)
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

float sumMatrix(
    int m,
    int n,
    float **in)
    {
    int i ,j;
    float sum=0;
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
    float **in)
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
    float **in)
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
    float ***u,
    float ** phi)
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
    float *****w,
    int iw,
    int jw,
    int kw,
    float **out)
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
    int M,
    int N,
    float ***u,
    int h,
    int p_s,
    float **kernel,
    int P_S,
    float **kernelk,
    int t_r,
    int s_r,
    int p_r,
    int k_r,
    int sw,
    float **phi,
    int s_s,
    float *****w)
    {

    int i, j, k, r, s, i0, j0, ii, jj, iii;

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<c;k++)
            {
                i0=i+t_r;
                j0=j+t_r;

                subMatrix   (p_r,M,N,phi,i0,j0,sphi);
                subMatrix   (k_r,M,N,phi,i0,j0,sphik);

                if(any(2*p_r+1,2*p_r+1,sphi) == 1)
                {

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
                                subMatrix3D(p_r,M,N,c,k,u,r,s,W2);
                                subMatrix3D(k_r,M,N,c,k,u,r,s,W2k);

                                subt2Matrix (2*p_r+1,2*p_r+1,W1,  W2  ,diff);
                                subt2Matrix (2*k_r+1,2*k_r+1,W1k, W2k ,diffk);

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
    int M,
    int N,
    float ***u,
    int h,
    int p_s,
    float **kernel,
    int P_S,
    float **kernelk,
    int t_r,
    int s_r,
    int p_r,
    int k_r,
    int sw,
    float **phi,
    float **PHI,
    int s_s,
    float *****w)
    {

    int i, j, k, r, s, i0, j0, ii, jj, iii, jjj;

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<c;k++)
            {
                i0=i+t_r;
                j0=j+t_r;

                subMatrix   (p_r,M,N,phi,i0,j0,sphi);
                subMatrix   (k_r,M,N,phi,i0,j0,sphik);

                if (PHI[i0][j0]==0)
                {
                    if(any(2*p_r+1,2*p_r+1,sphi) == 1)
                    {

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
                                    subMatrix3D(p_r,M,N,c,k,u,r,s,W2);
                                    subMatrix3D(k_r,M,N,c,k,u,r,s,W2k);

                                    subt2Matrix (2*p_r+1,2*p_r+1,W1,  W2  ,diff);
                                    subt2Matrix (2*k_r+1,2*k_r+1,W1k, W2k ,diffk);

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
    float ***u0,
    int M,
    int N,
    int h,
    int p_s,
    float **kernel,
    int P_S,
    float **kernelk,
    int t_r,
    int s_r,
    int p_r,
    int sw,
    float **phi,
    float **PHI,
    int s_s,
    float *****w,
    float lamda,
    float ***f0)
    {
    int i,j,k,i0,j0,iii,ii,jj;

    int last_count, step, done;

    int k_r=sw*p_r;

    initU(m,n,c,t_r,u0,u);
    initVars(p_s,P_S,s_s,M,N,m,n,c);

    last_count = sumMatrix(M,N,phi);
    printfFnc("INPAINTING \n");
    step =  1;
    done = -1;
    while(step < 99999 && done != 1)
    {
        padarray3d(m,n,c,t_r,u0,u);
        if(step==1)
        {
            printfFnc("Pierwsze obliczenie wagi: ");
            updateWeight(m,n,c,M,N,u,h,p_s,kernel,
                P_S,kernelk,t_r,s_r,p_r,k_r,sw,phi,s_s,w);
            printfFnc("OK. \n");
        }
        if(step>9 && step%10==0)
        {
            printfFnc("Step: %d \n", step);
            updatePhi(M,N,c,u,phi);
            updateWeight2(m,n,c,M,N,u,h,p_s,kernel,
                P_S,kernelk,t_r,s_r,p_r,k_r,sw,phi,PHI,s_s,w);
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

                        subMatrix(p_r,M,N,phi,i0,j0,sphi);
                        if(any(2*p_r+1,2*p_r+1,sphi)==1)

                        {

                            subWeight(m,n,c,s_s,w,i,j,k,subw);
                            subMatrix3D(s_r,M,N,c,k,u,i0,j0,subu);

                            prodMatrix(s_s,s_s,subw,subu,prod);
                            float sum_yw, sum_w;
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

        if (step>9 && step%10==0 && allOne(M,N,phi)==1)
        {
            printfFnc("DONE");
            done = 1;
        } else {
            last_count = sumMatrix(M,N,phi);
            step++;
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
    int P_S          = mxGetScalar(pmxIn[9]);
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

    int i,j,i0,j0,it,in4d,in2d,k,l,c_it,si,sj;

    
    initVars(p_s,P_S,s_s,M,N,m,n,c);

    c_it = 0;
    it   = 0;
    while(c_it<c)
    {
        while(it<(c_it+1)*m*n)
        {
            in2d=0;
            while(in2d<m*n)
            {
                u0[in2d%m][in2d/m][c_it] = (float)u0u[it];
                f0[in2d%m][in2d/m][c_it] = (float)f0u[it++];
                in2d++;
            }
        }
        c_it++;
    }

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            i0=i+t_r;
            j0=j+t_r;
            for(k=0;k<c;k++) {
                u[i0][j0][k]=u0[i][j][k];
            }
        }
    }

    for(i=0;i<p_s*p_s;i++)
    {
        kernel[i%p_s][i/p_s] = (float)kernelu[i];
    }

    for(i=0;i<P_S*P_S;i++)
    {
        kernelk[i%P_S][i/P_S] = (float)kernelku[i];
    }

    for(i=0;i<M*N;i++)
    {
        phi[i%M][i/M] = (float)phiu[i];
        PHI[i%M][i/M] = (float)PHIu[i];
    }

    solveNLCTV(m,n,c,u0,M,N,h,p_s,kernel,P_S,kernelk,t_r,
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

    printfFnc("Czyszczenie zmiennych: ");
    clearVars(p_s,P_S,s_s,M,N,m,n,c);
    printfFnc("OK. \n");
}
