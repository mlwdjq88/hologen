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
double * ExportPolygon(double * B, long rowLen, long  &num){
    double * Bs;
    num=ceil(((double)rowLen-2)/2);
//     std::cout << std::fixed <<num<<  '\n';
    if (rowLen<3){
        return 0;
    }
    else if (fmod(rowLen-1,2)==0){
        Bs= (double*) malloc(10*num * sizeof(double));
        Bs[0]=B[0];
        Bs[1]=B[rowLen];
        Bs[2]=B[0];
        Bs[3]=B[rowLen];
        Bs[4]=B[rowLen-1];
        Bs[5]=B[rowLen-1+rowLen];
        Bs[6]=B[1];
        Bs[7]=B[1+rowLen];
        Bs[8]=B[0];
        Bs[9]=B[rowLen];
        for( int i=1; i<num; i++){
            Bs[10*i]=B[i];
            Bs[10*i+1]=B[i+rowLen];
            Bs[10*i+2]=B[rowLen-i];
            Bs[10*i+3]=B[rowLen-i+rowLen];
            Bs[10*i+4]=B[rowLen-1-i];
            Bs[10*i+5]=B[rowLen-1-i+rowLen];
            Bs[10*i+6]=B[i+1];
            Bs[10*i+7]=B[i+1+rowLen];
            Bs[10*i+8]=B[i];
            Bs[10*i+9]=B[i+rowLen];
        }
    }
    else{
        Bs= (double*) malloc(10*num * sizeof(double));
        for( int i=0; i<num; i++){
            Bs[10*i]=B[i];
            Bs[10*i+1]=B[i+rowLen];
            Bs[10*i+2]=B[rowLen-1-i];
            Bs[10*i+3]=B[rowLen-1-i+rowLen];
            Bs[10*i+4]=B[rowLen-2-i];
            Bs[10*i+5]=B[rowLen-2-i+rowLen];
            Bs[10*i+6]=B[i+1];
            Bs[10*i+7]=B[i+1+rowLen];
            Bs[10*i+8]=B[i];
            Bs[10*i+9]=B[i+rowLen];
        }
    }
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
         free (C);
    }
    fclose(outputFile);
}
