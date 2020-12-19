#include <fstream>
#include <stdio.h>
#include <mex.h>
#include <math.h>
#include <iostream>
#include <cstring>
#define PI 3.14159265

double calDistance(double X1,double Y1,double X2,double Y2,double X0,double Y0){
    double A,B,C,d;
    A=Y2-Y1;
    B=X1-X2;
    C=X2*Y1-X1*Y2;
    d=fabs(A*X0+B*Y0+C)/sqrt(pow(A,2)+pow(B,2));
    return d;
}

double * downSampling(double * B, double lambda, double delta, double f,double * xs,double * ys,long & rowLen,long xsLen,double AccuCtrl){
    double * Bs,X1,X2,X0,Y1,Y2,Y0,d;
    double Rs,K,lv,localT;
    int k=0,t=0;
    long num=rowLen;
    for (int i=0;i<rowLen;i++){
        Rs=pow(xs[(long)((B[i+rowLen]-1)*xsLen+B[i]-1)],2)+pow(ys[(long)((B[i+rowLen]-1)*xsLen+B[i]-1)],2);
        K=K+sqrt(Rs+pow(f,2)+pow(delta/2,2))/lambda-PI/4;
    }
    K=K/rowLen;
    lv=pow(floor(K)+PI/4,2)*pow(lambda,2)-pow(f,2)-pow(delta/2,2);
    if (lv>=0){
        localT=sqrt(pow(ceil(K)+PI/4,2)*pow(lambda,2)-pow(f,2)-pow(delta/2,2))-sqrt(lv);
    }
    else{
        localT=sqrt(pow(ceil(K)+PI/4,2)*pow(lambda,2)-pow(f,2)-pow(delta/2,2));
    }
    for (int q=1;q<rowLen-1;q++){
        X1=xs[(long)((B[k+rowLen]-1)*xsLen+B[k]-1)];
        Y1=ys[(long)((B[k+rowLen]-1)*xsLen+B[k]-1)];
        X2=xs[(long)((B[q+1+rowLen]-1)*xsLen+B[q+1]-1)];
        Y2=ys[(long)((B[q+1+rowLen]-1)*xsLen+B[q+1]-1)];
        X0=xs[(long)((B[q+rowLen]-1)*xsLen+B[q]-1)];
        Y0=ys[(long)((B[q+rowLen]-1)*xsLen+B[q]-1)];
        d=calDistance(X1,Y1,X2,Y2,X0,Y0);
        if (d<=localT*AccuCtrl){
            B[q]=0;
            num--;
        }
        else{
            k=q;
        }
    }
    Bs= (double*) malloc(2*num * sizeof(double));
    
    for (int p=0;p<rowLen;p++){
        if (B[p]!=0){
            Bs[t]=xs[(long)((B[p+rowLen]-1)*xsLen+B[p]-1)];
            Bs[t+num]=ys[(long)((B[p+rowLen]-1)*xsLen+B[p]-1)];
            t++;
        }
    }
    rowLen=num;
    return Bs;
}



// Gateway function for MATLAB interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *C,*pB,*pC;
    mxArray *B,*Bs;
    long cellLen,rowLen,xsLen;
    cellLen          = mxGetM(prhs[0]);
//      cellLen=4;
    double lambda       = mxGetScalar(prhs[1]);
    double delta       = mxGetScalar(prhs[2]);
    double f           = mxGetScalar(prhs[3]);
    double *xs       = mxGetPr(prhs[4]);
    double *ys       = mxGetPr(prhs[5]);
    double AccuCtrl       = mxGetScalar(prhs[6]);
    xsLen          = mxGetM(prhs[4]);
    plhs[0]=mxCreateCellMatrix(cellLen, 1);
    for(int i=0; i<cellLen; i++){
        B=mxGetCell(prhs[0], i);
        rowLen          = mxGetM(B);
        pB = mxGetPr(B);
        C=downSampling(pB,lambda,delta,f,xs,ys,rowLen,xsLen,AccuCtrl);
        Bs=mxCreateDoubleMatrix(rowLen, 2,mxREAL);
        pC = mxGetPr(Bs);
        for(int j=0; j<rowLen; j++){
            pC[j]=C[j];
            pC[j+rowLen]=C[j+rowLen];
        }
        mxSetCell(plhs[0], i,Bs);
//         std::cout << std::fixed << i<<  '\n';
    }
}
