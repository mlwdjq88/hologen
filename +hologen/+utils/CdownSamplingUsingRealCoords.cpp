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

double * downSampling(double * B, double lambda, double delta, double f,long & rowLen,double AccuCtrl){
    double * Bs,X1,X2,X0,Y1,Y2,Y0,d;
    double Rs,K=0,lv,localT;
    int k=0,t=0;
    long num=rowLen;
    for (int i=0;i<rowLen;i++){
        Rs=pow(B[i],2)+pow(B[i+rowLen],2);
        K=K+(sqrt(Rs+pow(f,2))-f)/lambda-2;
    }
    K=K/rowLen;
    lv=pow(floor(K+2)*lambda+f,2)-pow(f,2);
    if (lv>=0){
        localT=sqrt(pow(ceil(K+2)*lambda+f,2)-pow(f,2))-sqrt(lv);
    }
    else{
        localT=sqrt(pow(ceil(K+2)*lambda+f,2)-pow(f,2));
    }
//      std::cout << std::fixed <<K<<'\t'<<localT<<  '\n';
    for (int q=1;q<rowLen-1;q++){
        X1=B[k];
        Y1=B[k+rowLen];
        X2=B[q+1];
        Y2=B[q+1+rowLen];
        X0=B[q];
        Y0=B[q+rowLen];
        d=calDistance(X1,Y1,X2,Y2,X0,Y0);
//                  
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
            Bs[t]=B[p];
            Bs[t+num]=B[p+rowLen];
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
    long cellLen,rowLen;
    cellLen          = mxGetM(prhs[0]);
//      cellLen=4;
    double lambda       = mxGetScalar(prhs[1]);
    double delta       = mxGetScalar(prhs[2]);
    double f           = mxGetScalar(prhs[3]);
    double AccuCtrl       = mxGetScalar(prhs[4]);
    plhs[0]=mxCreateCellMatrix(cellLen, 1);
    for(int i=0; i< cellLen ; i++){
        B=mxGetCell(prhs[0], i);
        rowLen          = mxGetM(B);
        pB = mxGetPr(B);
        C=downSampling(pB,lambda,delta,f,rowLen,AccuCtrl);
        Bs=mxCreateDoubleMatrix(rowLen, 2,mxREAL);
        pC = mxGetPr(Bs);
        for(int j=0; j<rowLen; j++){
            pC[j]=C[j];
            pC[j+rowLen]=C[j+rowLen];
        }
        mxSetCell(plhs[0], i,Bs);
        free(C);
//         std::cout << std::fixed << i<<  '\n';
    }
}
