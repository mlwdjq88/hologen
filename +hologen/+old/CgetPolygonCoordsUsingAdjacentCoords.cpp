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
double * ExportPolygon(double * B, long rowLen, long  &num, double offsetx,double offsety,double R,double RStop){
    double * Bs;
    num=ceil(((double)rowLen-2)/2);
     
    if (rowLen<3){
        return 0;
    }
    else if ( pow(B[0],2)+pow(B[rowLen],2)>pow(R,2)||pow(B[0],2)+pow(B[rowLen],2)<pow(RStop,2)){
        num=0;
        return 0;
    }
    else if (fmod(rowLen-1,2)==0){
        Bs= (double*) malloc(10*num * sizeof(double));
        Bs[0]=B[0]+offsetx;
        Bs[1]=B[rowLen]+offsety;
        Bs[2]=B[0]+offsetx;
        Bs[3]=B[rowLen]+offsety;
        Bs[4]=B[rowLen-1]+offsetx;
        Bs[5]=B[rowLen-1+rowLen]+offsety;
        Bs[6]=B[1]+offsetx;
        Bs[7]=B[1+rowLen]+offsety;
        Bs[8]=B[0]+offsetx;
        Bs[9]=B[rowLen]+offsety;
        for( int i=1; i<num; i++){
            Bs[10*i]=B[i]+offsetx;
            Bs[10*i+1]=B[i+rowLen]+offsety;
            Bs[10*i+2]=B[rowLen-i]+offsetx;
            Bs[10*i+3]=B[rowLen-i+rowLen]+offsety;
            Bs[10*i+4]=B[rowLen-1-i]+offsetx;
            Bs[10*i+5]=B[rowLen-1-i+rowLen]+offsety;
            Bs[10*i+6]=B[i+1]+offsetx;
            Bs[10*i+7]=B[i+1+rowLen]+offsety;
            Bs[10*i+8]=B[i]+offsetx;
            Bs[10*i+9]=B[i+rowLen]+offsety;
        }
    }
    else{
        Bs= (double*) malloc(10*num * sizeof(double));
        for( int i=0; i<num; i++){
            Bs[10*i]=B[i]+offsetx;
            Bs[10*i+1]=B[i+rowLen]+offsety;
            Bs[10*i+2]=B[rowLen-1-i]+offsetx;
            Bs[10*i+3]=B[rowLen-1-i+rowLen]+offsety;
            Bs[10*i+4]=B[rowLen-2-i]+offsetx;
            Bs[10*i+5]=B[rowLen-2-i+rowLen]+offsety;
            Bs[10*i+6]=B[i+1]+offsetx;
            Bs[10*i+7]=B[i+1+rowLen]+offsety;
            Bs[10*i+8]=B[i]+offsetx;
            Bs[10*i+9]=B[i+rowLen]+offsety;
        }
    }
 //   std::cout << std::fixed <<offsetx<<  '\n';
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
    double offsetx       = mxGetScalar(prhs[7]);
    double offsety       = mxGetScalar(prhs[8]);
    double R       = mxGetScalar(prhs[9]);
    double RStop       = mxGetScalar(prhs[10]);
    long num;
    char filename[20]="test";
    strcat(strcat(strcat(strcat(filename,num2str(p)),"&"),num2str(q)),".cor");
//     strcat(filename,num2str(p));
    //  std::cout << std::fixed <<offset<<  '\n';
    char *cCoords;
    if ((outputFile = fopen(filename, "wb")) == NULL)
        printf("Cannot open file.\n");
    for(int i=0; i<cellLen; i++){
        num=0;
        B=mxGetCell(prhs[0], i);
        rowLen          = mxGetM(B);
        pB = mxGetPr(B);
        C=ExportPolygon(pB,(long)rowLen,num,offsetx,offsety,R,RStop);
        if (num>0){
            E=new long [16*num];
                for (int j=0; j<num; j++){
                    E[16*j]=polyPre[0];
                    E[16*j+1]=polyPre[1];
                    E[16*j+2]=polyPre[2];
                    E[16*j+3]=polyPre[3];
                    E[16*j+4]=polyForm;
                    for (int k=0; k<10; k++){
                        E[16*j+5+k]=C[10*j+k];
//                          std::cout << std::fixed << C[10*j+2*k]<<  '\n';
                    }
                    E[16*j+15]=polyPost;
//                     fwrite(E, sizeof(long), 4*16*num/5, outputFile);
//                       std::cout << std::fixed <<  E[16*j+15+16*5*m]<<'\t' <<16*j+15+16*5*m<<'\n';
                }
//             fwrite(E, sizeof(long), 4*16*num, outputFile);
             cCoords=new char [4*16*num];
             encodePoly32(E, cCoords,16*num);
             fwrite(cCoords, sizeof(char), 4*16*num, outputFile);
             delete cCoords;
             delete E;
//                std::cout << std::fixed << atan2(-1,0)<<  '\n';
        }
         free (C);
    }
    fclose(outputFile);
}
