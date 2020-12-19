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

double dmax(double * A, int Ibegin, int Iend){
    double temp;
    int i;
    for(i=Ibegin; i<Iend; i++){
        if (A[i]>A[i+1]){
            temp=A[i];
            A[i]=A[i+1];
            A[i+1]=temp;
        }
    }
    return A[i];
}

double dmin(double * A, int Ibegin, int Iend){
    double temp;
    int i;
    for(i=Ibegin; i<Iend; i++){
        if (A[i]<A[i+1]){
            temp=A[i];
            A[i]=A[i+1];
            A[i+1]=temp;
        }
    }
    return A[i];
}

double * ExportPolygon(double * B, int dire, long rowLen, long  &num){
    double * Bs;
    int i;
    int k=0;
    int t=1;
    double *D= new double [10*rowLen];
    sort(B,dire,rowLen);
//     for( i=(dire-1)*rowLen+1; i<dire*rowLen-1; i++){
//         if (B[i]<B[i+1]){
//             k++;
//         }
//     }
    if (dire==1){
        for( i=(dire-1)*rowLen; i<dire*rowLen-1; i++){
            if (B[i]<B[i+1]){
                if (k==0){
                    D[2*num]=B[i];
                    D[2*num+1]=B[i+rowLen];
                    num++;
                    t++;
                }
                else{
                    if (t==3&&D[2*num-3]>D[2*num-1]){
                        D[2*num]=B[i];
                        D[2*num+1]=dmin(B,i-k+rowLen,i+rowLen);
                        D[2*num+2]=B[i];
                        D[2*num+3]=dmax(B,i-k+rowLen,i+rowLen);
                        
                    }
                    else{
                        D[2*num]=B[i];
                        D[2*num+1]=dmax(B,i-k+rowLen,i+rowLen);
                        D[2*num+2]=B[i];
                        D[2*num+3]=dmin(B,i-k+rowLen,i+rowLen);
                        
                    }
                    
                    t=t+2;
                    num=num+2;
                }
                k=0;
            }
            else if(i==dire*rowLen-2){
                if (t==3&&D[2*num-3]>D[2*num-1]){
                    D[2*num]=B[i];
                    D[2*num+1]=dmin(B,i-k+rowLen,i+rowLen+1);
                    D[2*num+2]=B[i];
                    D[2*num+3]=dmax(B,i-k+rowLen,i+rowLen+1);
                }
                else{
                    D[2*num]=B[i];
                    D[2*num+1]=dmax(B,i-k+rowLen,i+rowLen+1);
                    D[2*num+2]=B[i];
                    D[2*num+3]=dmin(B,i-k+rowLen,i+rowLen+1);
                }
                t=t+2;
                num=num+2;
            }
            else{
                k++;
            }
            if (t==4){
                D[2*num]=D[2*num-6];
                D[2*num+1]=D[2*num-5];
                D[2*num+2]=D[2*num-6];
                D[2*num+3]=D[2*num-5];
                num=num+2;
                t=3;
                if (i!=dire*rowLen-2){
                    D[2*num]=D[2*num-8];
                    D[2*num+1]=D[2*num-7];
                    D[2*num+2]=D[2*num-6];
                    D[2*num+3]=D[2*num-5];
                    num=num+2;
                }
            }
            else if (t==5){
                D[2*num]=D[2*num-8];
                D[2*num+1]=D[2*num-7];
                num=num+1;
                t=3;
                if (i!=dire*rowLen-2){
                    D[2*num]=D[2*num-6];
                    D[2*num+1]=D[2*num-5];
                    D[2*num+2]=D[2*num-4];
                    D[2*num+3]=D[2*num-3];
                    num=num+2;
                }
            }
        }
    }
    else{
        for( i=(dire-1)*rowLen; i<dire*rowLen-1; i++){
            if (B[i]<B[i+1]){
                if (k==0){
                    D[2*num]=B[i-rowLen];
                    D[2*num+1]=B[i];
                    num++;
                    t++;
                }
                else{
                    if (t==3&&D[2*num-4]>D[2*num-2]){
                        D[2*num]=dmin(B,i-k-rowLen,i-rowLen);
                        D[2*num+1]=B[i];
                        D[2*num+2]=dmax(B,i-k-rowLen,i-rowLen);
                        D[2*num+3]=B[i];
                    }
                    else{
                        D[2*num]=dmax(B,i-k-rowLen,i-rowLen);
                        D[2*num+1]=B[i];
                        D[2*num+2]=dmin(B,i-k-rowLen,i-rowLen);
                        D[2*num+3]=B[i];
                    }
                    t=t+2;
                    num=num+2;
                }
                k=0;
            }
            else if(i==dire*rowLen-2){
                if (t==3&&D[2*num-4]>D[2*num-2]){
                    D[2*num]=dmin(B,i-k-rowLen,i-rowLen+1);
                    D[2*num+1]=B[i];
                    D[2*num+2]=dmax(B,i-k-rowLen,i-rowLen+1);
                    D[2*num+3]=B[i];
                }
                else{
                    D[2*num]=dmax(B,i-k-rowLen,i-rowLen+1);
                    D[2*num+1]=B[i];
                    D[2*num+2]=dmin(B,i-k-rowLen,i-rowLen+1);
                    D[2*num+3]=B[i];
                }
                t=t+2;
                num=num+2;
            }
            else{
                k++;
            }
            if (t==4){
                D[2*num]=D[2*num-6];
                D[2*num+1]=D[2*num-5];
                D[2*num+2]=D[2*num-6];
                D[2*num+3]=D[2*num-5];
                num=num+2;
                t=3;
                if (i!=dire*rowLen-2){
                    D[2*num]=D[2*num-8];
                    D[2*num+1]=D[2*num-7];
                    D[2*num+2]=D[2*num-6];
                    D[2*num+3]=D[2*num-5];
                    num=num+2;
                }
            }
            else if (t==5){
                D[2*num]=D[2*num-8];
                D[2*num+1]=D[2*num-7];
                num=num+1;
                t=3;
                if (i!=dire*rowLen-2){
                    D[2*num]=D[2*num-6];
                    D[2*num+1]=D[2*num-5];
                    D[2*num+2]=D[2*num-4];
                    D[2*num+3]=D[2*num-3];
                    num=num+2;
                }
            }
        }
    }
    Bs= (double*) malloc(2*num * sizeof(double));
    for( i=0; i<num; i++){
        Bs[2*i]=D[2*i];
        Bs[2*i+1]=D[2*i+1];
    }
    delete D;
//         for( i=0; i<num; i++){
//         std::cout << std::fixed << Bs[2*i] <<'\t' <<Bs[2*i+1]<<'\n';
//     }
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
//     unsigned char gdsPost[8];
//     unsigned char polyPre[16];
//     unsigned char polyPost[4];
//     unsigned char polyForm[4];
//     int layerNumber =1;
//     char * fileName ="test";
    double *C,*pB;
    long *E;
    mxArray *B;
    mwSize cellLen,rowLen;
    cellLen          = mxGetM(prhs[0]);
//      cellLen=4;
    int dire       = mxGetScalar(prhs[1]);
    double *polyPre       = mxGetPr(prhs[2]);
    double polyForm       = mxGetScalar(prhs[3]);
    double polyPost       = mxGetScalar(prhs[4]);
    double *xs       = mxGetPr(prhs[5]);
    double *ys       = mxGetPr(prhs[6]);
    long Ny       = mxGetScalar(prhs[7]);
    int p       = mxGetScalar(prhs[8]);
    int q       = mxGetScalar(prhs[9]);
    long len;
    long num;
    char filename[20]="test";
    strcat(strcat(strcat(strcat(filename,num2str(p)),"&"),num2str(q)),".cor");
//     strcat(filename,num2str(p));
//      std::cout << std::fixed <<filename<<  '\n';
//     char *cCoords;
    int flip[8]={1,1,1,-1,-1,-1,-1,1};
    if ((outputFile = fopen(filename, "wb")) == NULL)
                printf("Cannot open file.\n");
    for(int i=0; i<cellLen; i++){
        num=0;
        B=mxGetCell(prhs[0], i);
        rowLen          = mxGetM(B);
        pB = mxGetPr(B);
        C=ExportPolygon(pB,dire,(long)rowLen,num);
        if (num>=5){
            E=new long [4*16*num/5];
            for(int m=0; m<4; m++){
                for (int j=0; j<num/5; j++){
                    E[16*j+16*num/5*m]=polyPre[0];
                    E[16*j+1+16*num/5*m]=polyPre[1];
                    E[16*j+2+16*num/5*m]=polyPre[2];
                    E[16*j+3+16*num/5*m]=polyPre[3];
                    E[16*j+4+16*num/5*m]=polyForm;
                    for (int k=0; k<5; k++){
                        E[16*j+5+2*k+16*num/5*m]=flip[2*m]*xs[(long)((2*Ny+1)*(C[10*j+2*k+1]-1)+C[10*j+2*k]-1)];
                        E[16*j+6+2*k+16*num/5*m]=flip[2*m+1]*ys[(long)((2*Ny+1)*(C[10*j+2*k+1]-1)+C[10*j+2*k]-1)];
                    }
                    E[16*j+15+16*num/5*m]=polyPost;
//                     fwrite(E, sizeof(long), 4*16*num/5, outputFile);
//                       std::cout << std::fixed <<  E[16*j+15+16*5*m]<<'\t' <<16*j+15+16*5*m<<'\n';
                }
            }
//             cCoords=new char [4*4*16*num/5];
//             encodePoly32(E, cCoords,4*16*num/5);
            fwrite(E, sizeof(long), 4*16*num/5, outputFile);
//             delete cCoords;
            delete E;
//               std::cout << std::fixed << "cellLen " << num<< " rowLen "  << 4*16*num/5<<  '\n';
        }
    }
    fclose(outputFile);
    
    
}
