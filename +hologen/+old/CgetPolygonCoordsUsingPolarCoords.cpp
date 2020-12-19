/* InitGDS.cpp, Wenhua, 10/28/2018
 * this function is used to initialize GDS file
 *
 */

#include <fstream>
#include <stdio.h>
#include <mex.h>
#include <math.h>
#include <iostream>
#include <cstring>
// #include <sstream>
// #include <string>

char * num2str(int i)
{
    char tmp[100];
    sprintf(tmp, "%d", i);
    return tmp;
}

void avg(double *Pb,int n,double *pm)
{
    pm[0]=0;
    pm[1]=0;  
    for(int i=0; i<n; i++){
        pm[0]=pm[0]+Pb[i]/(double)n;
        pm[1]=pm[1]+Pb[i+n]/(double)n;
    }
}

void sort(double * A, int dire, long rowLen){
//     double *B;
    double temp;
//     B=new double [Iend-Ibegin+1];
    int i,j;
    if( dire==1){
        for(i=0; i<rowLen-1; i++){
            for( j=0; j<rowLen-1-i; j++){
                if (A[j]>A[j+1]){
                    temp=A[j];
                    A[j]=A[j+1];
                    A[j+1]=temp;
                    temp=A[j+rowLen];
                    A[j+rowLen]=A[j+1+rowLen];
                    A[j+1+rowLen]=temp;
                }
            }
        }
    }
    else{
        for(i=rowLen; i<2*rowLen-1; i++){
            for( j=rowLen; j<3*rowLen-1-i; j++){
                if (A[j]>A[j+1]){
                    temp=A[j];
                    A[j]=A[j+1];
                    A[j+1]=temp;
                    temp=A[j-rowLen];
                    A[j-rowLen]=A[j+1-rowLen];
                    A[j+1-rowLen]=temp;
                }
            }
        }
    }
}


double * ExportPolygon(double * B, long rowLen, long  &num){
    double * Bs,*P,*Ps,*BPs,*BPsm;
    P= (double*) malloc(2*rowLen * sizeof(double));
    for( int i=0; i<rowLen; i++){
        P[i]=atan2(B[i+rowLen],B[i]);
        P[i+rowLen]=sqrt(pow(B[i+rowLen],2)+pow(B[i],2));
    }
    sort(P,1,rowLen);
    num=ceil((rowLen-2)/2);
    if (rowLen<4){
        return 0;
    }
    else if (fmod(rowLen,2)==0){
        Bs= (double*) malloc(10*num * sizeof(double));
        BPs= (double*) malloc(8 * sizeof(double));
        BPsm= (double*) malloc(2 * sizeof(double));
        Ps= (double*) malloc(8 * sizeof(double));
        for( int i=0; i<num; i++){
            for( int j=0; j<4; j++){
                Ps[j]=P[2*i+j];
                Ps[j+4]=P[2*i+j+rowLen];
                BPs[j]=Ps[j+4]*cos(Ps[j]);
                BPs[j+4]=Ps[j+4]*sin(Ps[j]);
            }
            
            avg(BPs,4,BPsm);
            for( int j=0; j<4; j++){
                Ps[j]=atan2(BPs[j+4]-BPsm[1],BPs[j]-BPsm[0]);
                Ps[j+4]=sqrt(pow(BPs[j+4]-BPsm[1],2)+pow(BPs[j]-BPsm[0],2));
            }
            sort(Ps,1,4);
            for( int j=0; j<4; j++){
                Bs[10*i+2*j]=Ps[j+4]*cos(Ps[j])+BPsm[0];
                Bs[10*i+2*j+1]=Ps[j+4]*sin(Ps[j])+BPsm[1];
            }
            Bs[10*i+8]=Ps[4]*cos(Ps[0])+BPsm[0];
            Bs[10*i+9]=Ps[4]*sin(Ps[0])+BPsm[1];
        }
        free (BPs);
        free (BPsm);
        free (Ps);
    }
    else{
        Bs= (double*) malloc(10*num * sizeof(double));
        BPs= (double*) malloc(8 * sizeof(double));
        BPsm= (double*) malloc(2 * sizeof(double));
        Ps= (double*) malloc(8 * sizeof(double));
        for( int j=0; j<3; j++){
            Ps[j]=P[j];
            Ps[j+3]=P[j+rowLen];
            BPs[j]=Ps[j+3]*cos(Ps[j]);
            BPs[j+3]=Ps[j+3]*sin(Ps[j]);
        }
        for( int j=0; j<3; j++){
            Bs[2*j]=BPs[j];
            Bs[2*j+1]=BPs[j+3];
        }
        Bs[6]=BPs[0];
        Bs[7]=BPs[3];
        Bs[8]=BPs[0];
        Bs[9]=BPs[3];
        for( int i=1; i<num; i++){
            for( int j=0; j<4; j++){
                Ps[j]=P[2*i+j];
                Ps[j+4]=P[2*i+j+rowLen];
                BPs[j]=Ps[j+4]*cos(Ps[j]);
                BPs[j+4]=Ps[j+4]*sin(Ps[j]);
            }
                            
            avg(BPs,4,BPsm);
//             BPsm[0]=0;
//             BPsm[1]=0;
//             std::cout << std::fixed<<num<<'\t'<<BPsm[0]<< '\n';
            for( int j=0; j<4; j++){
                Ps[j]=atan2(BPs[j+4]-BPsm[1],BPs[j]-BPsm[0]);
                Ps[j+4]=sqrt(pow(BPs[j+4]-BPsm[1],2)+pow(BPs[j]-BPsm[0],2));
//                 std::cout << std::fixed<<BPsm[0]<< '\n';
            }
            sort(Ps,1,4);
            for( int j=0; j<4; j++){
                Bs[10*i+2*j]=Ps[j+4]*cos(Ps[j])+BPsm[0];
                Bs[10*i+2*j+1]=Ps[j+4]*sin(Ps[j])+BPsm[1];
            }
            Bs[10*i+8]=Ps[4]*cos(Ps[0])+BPsm[0];
            Bs[10*i+9]=Ps[4]*sin(Ps[0])+BPsm[1];
            
        }
        free (BPs);
        free (Ps);
        free (BPsm);
    }
    free (P);
    return Bs;
}
void encode32(long aCoord, char * cPart){
    cPart[0] = (aCoord >> 24) & 255;
    cPart[1] = (aCoord >> 16) & 255;
    cPart[2] = (aCoord >> 8) & 255;
    cPart[3] = (aCoord) & 255;
}

void encodePoly32(long * coords, char * cCoords,long num){
    char cPart[4];
    for (int k = 0; k < num; k++){
        encode32(coords[k], cPart);
        for (int m = 0; m < 4; m++){
            cCoords[k*4 + m] = cPart[m];
        }
    }
}

// Gateway function for MATLAB interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    FILE * outputFile = NULL;
    double *C,*pB;
    long *E;
    mxArray *B;
    mwSize cellLen,rowLen;
    cellLen          = mxGetM(prhs[0]);
    double *polyPre       = mxGetPr(prhs[1]);
    double polyForm       = mxGetScalar(prhs[2]);
    double polyPost       = mxGetScalar(prhs[3]);
    long Ny       = mxGetScalar(prhs[4]);
    int p       = mxGetScalar(prhs[5]);
    int q       = mxGetScalar(prhs[6]);
    long num;
    char filename[20]="test";
    strcat(strcat(strcat(strcat(filename,num2str(p)),"&"),num2str(q)),".cor");
//     strcat(filename,num2str(p));
//      std::cout << std::fixed <<filename<<  '\n';
    char *cCoords;
    int flip[8]={1,1,1,-1,-1,-1,-1,1};
    if ((outputFile = fopen(filename, "wb")) == NULL)
        printf("Cannot open file.\n");
    for(int i=0; i<cellLen; i++){
        num=0;
        B=mxGetCell(prhs[0], i);
        rowLen          = mxGetM(B);
        pB = mxGetPr(B);
        C=ExportPolygon(pB,(long)rowLen,num);
        if (num>0){
            E=new long [4*16*num];
            for(int m=0; m<4; m++){
                for (int j=0; j<num; j++){
                    E[16*j+16*num*m]=polyPre[0];
                    E[16*j+1+16*num*m]=polyPre[1];
                    E[16*j+2+16*num*m]=polyPre[2];
                    E[16*j+3+16*num*m]=polyPre[3];
                    E[16*j+4+16*num*m]=polyForm;
                    for (int k=0; k<5; k++){
                        E[16*j+5+2*k+16*num*m]=flip[2*m]*C[10*j+2*k];
                        E[16*j+6+2*k+16*num*m]=flip[2*m+1]*C[10*j+2*k+1];
//                          std::cout << std::fixed << C[10*j+2*k]<<  '\n';
                    }
                    E[16*j+15+16*num*m]=polyPost;
//                     fwrite(E, sizeof(long), 4*16*num/5, outputFile);
//                       std::cout << std::fixed <<  E[16*j+15+16*5*m]<<'\t' <<16*j+15+16*5*m<<'\n';
                }
            }
//             fwrite(E, sizeof(long), 4*16*num, outputFile);
            cCoords=new char [4*4*16*num];
            encodePoly32(E, cCoords,4*16*num);
            fwrite(cCoords, sizeof(char), 4*4*16*num, outputFile);
            delete cCoords;
            delete E;
//                std::cout << std::fixed << atan2(-1,0)<<  '\n';
        }
//         free (C);
    }
    fclose(outputFile);
}
