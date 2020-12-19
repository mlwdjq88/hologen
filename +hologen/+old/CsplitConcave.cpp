#include <fstream>
#include <stdio.h>
#include <mex.h>
#include <math.h>
#include <iostream>
#include <cstring>

double mean(double *A,long num){
    double s=0;
    for(long i=0;i<num;i++){
        s+=A[i];
    }
    s=s/num;
    return s;
}

long removeNULL(double *x,long num){
    double y[6000];
    long k=0;
    for(int i=0;i<num;i++){
        if (x[i]!=NULL){
            y[k]=x[i];
            k++;
        }
    }
    for(int i=0;i<k;i++){
            x[i]=y[i];
        }
    return k;
}

void circshift(double *x,int shift,int num,int N){
    double y[6000];
    for (int i=0;i<N;i++){
        for(int j=0;j<num;j++){
            if (j-shift>=0&&j-shift<num)
                y[j+i*num]=x[j-shift+i*num];
            else if(j-shift<0)
                y[j+i*num]=x[j-shift+num+i*num];
            else if(j-shift>=num)
                y[j+i*num]=x[j-shift-num+i*num];
        }
    }
    for(int i=0;i<N*num;i++)
        x[i]=y[i];
}

double calAngle(double x0,double y0,double x1,double y1,double x2,double y2){
    double angle,y,x;
    y=(x1-x0)*(y2-y1)-(x2-x1)*(y1-y0);
    x=(x1-x0)*(x2-x1)+(y1-y0)*(y2-y1);
    angle=atan2(y,x);
    return angle;
}


void splitConcave(double *matBs,int *pB,int cellLen,double *matCs,int *pC,int &nC){
    int rowLen0,rowLen,convex[3000],convexLen,flag;
    double angle[3000],sign,matB[6000];
    pC[0]=0;
    nC=0;
    for (int i=0;i<cellLen;i++){
        convexLen=0;
//          if (i==1){
//          std::cout << std::fixed <<i<<'\t'<<nC<<'\n';break;}
        rowLen=(pB[i+1]-pB[i])/2;
        if( rowLen>3000){
            std::cout << std::fixed <<"allocated array for split is too small"<<'\n';
            continue;
        }
        angle[0]=calAngle(matBs[pB[i]+rowLen-1],matBs[pB[i]+2*rowLen-1],matBs[pB[i]],matBs[pB[i]+rowLen],matBs[pB[i]+1],matBs[pB[i]+1+rowLen]);
        angle[rowLen-1]=calAngle(matBs[pB[i]+rowLen-2],matBs[pB[i]+2*rowLen-2],matBs[pB[i]+rowLen-1],matBs[pB[i]+2*rowLen-1],matBs[pB[i]],matBs[pB[i]+rowLen]);
        matB[0]=matBs[pB[i]];
        matB[rowLen]=matBs[pB[i]+rowLen];
        matB[rowLen-1]=matBs[pB[i]+rowLen-1];
        matB[2*rowLen-1]=matBs[pB[i]+2*rowLen-1];
        for (int j=1;j<rowLen-1;j++){
            angle[j]=calAngle(matBs[pB[i]+j-1],matBs[pB[i]+j-1+rowLen],matBs[pB[i]+j],matBs[pB[i]+j+rowLen],matBs[pB[i]+j+1],matBs[pB[i]+j+1+rowLen]);
            matB[j]=matBs[pB[i]+j];
            matB[j+rowLen]=matBs[pB[i]+j+rowLen];
        }
        if(mean(angle,rowLen)>0)
            sign=1;
        else
            sign=-1;
        for (int j=0;j<rowLen;j++){
            if (angle[j]*sign<0){
                convex[convexLen]=j;
                convexLen++;
            }
        }
        while (convexLen>0&&rowLen>2){
            circshift(matB,1-convex[0]-1,rowLen,2);
            circshift(angle,1-convex[0]-1,rowLen,1);
            if(angle[1]*sign>0){
                matCs[pC[nC]]=matB[0];
                matCs[pC[nC]+3]=matB[rowLen];
                matCs[pC[nC]+1]=matB[1];
                matCs[pC[nC]+4]=matB[rowLen+1];
                matCs[pC[nC]+2]=matB[2];
                matCs[pC[nC]+5]=matB[rowLen+2];
                pC[nC+1]=pC[nC]+6;
                nC++;
                matB[1]=NULL;
                matB[rowLen+1]=NULL;
                angle[1]=NULL;
                rowLen0=removeNULL(matB,2*rowLen)/2;
                rowLen=removeNULL(angle,rowLen);
                angle[0]=calAngle(matB[rowLen-1],matB[2*rowLen-1],matB[0],matB[rowLen],matB[1],matB[rowLen+1]);
                angle[1]=calAngle(matB[0],matB[rowLen],matB[1],matB[rowLen+1],matB[2],matB[rowLen+2]);
            }
            else if (angle[rowLen-1]*sign>0){
                matCs[pC[nC]]=matB[0];
                matCs[pC[nC]+3]=matB[rowLen];
                matCs[pC[nC]+1]=matB[rowLen-2];
                matCs[pC[nC]+4]=matB[2*rowLen-2];
                matCs[pC[nC]+2]=matB[rowLen-1];
                matCs[pC[nC]+5]=matB[2*rowLen-1];
                pC[nC+1]=pC[nC]+6;
                nC++;
                matB[rowLen-1]=NULL;
                matB[2*rowLen-1]=NULL;
                angle[rowLen-1]=NULL;
                rowLen0=removeNULL(matB,2*rowLen)/2;
                rowLen=removeNULL(angle,rowLen);
                angle[0]=calAngle(matB[rowLen-1],matB[2*rowLen-1],matB[0],matB[rowLen],matB[1],matB[rowLen+1]);
                angle[rowLen-1]=calAngle(matB[rowLen-2],matB[2*rowLen-2],matB[rowLen-1],matB[2*rowLen-1],matB[0],matB[rowLen]);
            }
            else{
                flag=0;
                for (int q=0;q<rowLen;q++)
                    if (angle[q]*sign<=0)
                        flag++;
                if (flag==rowLen)
                    sign=sign*(-1);
                circshift(matB,-1,rowLen,2);
            	circshift(angle,-1,rowLen,1);
            }
            convexLen=0;
            for (int j=0;j<rowLen;j++){
                if (angle[j]*sign<0){
                    convex[convexLen]=j;
                    convexLen++;
                }
            }
        }
        for (int j=0;j<rowLen;j++){
            matCs[pC[nC]+j]=matB[j];
            matCs[pC[nC]+j+rowLen]=matB[j+rowLen];
        }
        pC[nC+1]=pC[nC]+2*rowLen;
        nC++;
    }
}


// Gateway function for MATLAB interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int cellLen          = mxGetM(prhs[0]);
    mxArray *cellB,*Bs;
    double *matB,*matBs,*matCs,*matC;
    int totalLen=0,*pB,*pC,nC=0,rowLen;
    pB=(int*) malloc( (cellLen+1) * sizeof(int));
    pB[0]=0;
    for(int i=0; i< cellLen ; i++){
        cellB=mxGetCell(prhs[0], i);
        totalLen += mxGetM(cellB);
        pB[i+1] =2*totalLen;
    }
    matBs=(double*) malloc( 2*totalLen * sizeof(double));
    matCs=(double*) malloc( 3*2*totalLen * sizeof(double));
    pC=(int*) malloc( totalLen * sizeof(int));
    for(int i=0; i< cellLen ; i++){
        cellB=mxGetCell(prhs[0], i);
        rowLen= mxGetM(cellB);
        matB = mxGetPr(cellB);
        for(int j=0; j<rowLen; j++){
            matBs[pB[i]+j]=matB[j];
            matBs[pB[i]+j+rowLen]=matB[j+rowLen];
        }
    }
    splitConcave(matBs,pB,cellLen,matCs,pC,nC);
    plhs[0]=mxCreateCellMatrix(nC, 1);
    for(int i=0; i< nC ; i++){
        rowLen=(pC[i+1]-pC[i])/2;
        Bs=mxCreateDoubleMatrix(rowLen, 2,mxREAL);
        matC = mxGetPr(Bs);
        for(int j=0; j<rowLen; j++){
            matC[j]=matCs[pC[i]+j];
            matC[j+rowLen]=matCs[pC[i]+j+rowLen];
        }
        mxSetCell(plhs[0], i,Bs);
    }
    free(matBs);
    free(matCs);
    free(pB);
    free(pC);
}
