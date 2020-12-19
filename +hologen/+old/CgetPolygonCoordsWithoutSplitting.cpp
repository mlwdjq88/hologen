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
    if (rowLen<3){
        return 0;
    }
    else if ( pow(B[0],2)+pow(B[rowLen],2)>pow(R,2)||pow(B[0],2)+pow(B[rowLen],2)<pow(RStop,2)){
        num=0;
        return 0;
    }
    else if (B[0]!=B[rowLen-1]||B[rowLen]!=B[2*rowLen-1]){
        Bs= (double*) malloc((rowLen+2) * sizeof(double));
        num=rowLen+1;
        for( int i=0; i< rowLen; i++){
            Bs[i]=B[i]+offsetx;
            Bs[i+num]=B[i+rowLen]+offsety;
        }
         Bs[rowLen]=Bs[0];
         Bs[2*rowLen+1]=Bs[num];
    }
    else{
        Bs= (double*) malloc((rowLen) * sizeof(double));
        num=rowLen;
        for( int i=0; i< rowLen; i++){
            Bs[i]=B[i]+offsetx;
            Bs[i+num]=B[i+rowLen]+offsety;
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
            E=new long [2*num+6];
//                 for (int j=0; j<num; j++){
                    E[0]=polyPre[0];
                    E[1]=polyPre[1];
                    E[2]=polyPre[2];
                    E[3]=polyPre[3];
                    E[4]=polyForm;
                    for (int k=0; k<2*num; k++){
                        E[5+k]=C[k];
                    }
                    E[2*num+5]=polyPost;
//                     fwrite(E, sizeof(long), 4*16*num/5, outputFile);
//                       std::cout << std::fixed <<  E[16*j+15+16*5*m]<<'\t' <<16*j+15+16*5*m<<'\n';
//                 }
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
