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
#include <stdlib.h>
// #include <sstream>
// #include <string>

char * num2str(int i)
{
	char tmp[100];
	sprintf(tmp, "%d", i);
	return tmp;
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

void renderGDS(FILE * outputFile, unsigned char * gdsPost){
    fwrite(gdsPost, sizeof(char), 8, outputFile);
    fclose(outputFile);
}

void initGDS(FILE * outputFile, unsigned char * gdsPost, char * filename,int clen){
    int gdspreamble[102] = { 0, 6, 0, 2, 0, 7, 0, 28, 1, 2, 230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
        230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 10, 2, 6, 110, 111, 110, 97,
        109, 101, 0, 20, 3, 5, 61, 104, 219, 139, 172, 113, 12, 180, 56, 109,
        243, 127, 103, 94, 246, 236, 0, 28, 5, 2, 0, 114, 0, 4, 0, 17, 0, 13,
        0, 22, 0, 56, 0, 114, 0, 4, 0, 17, 0, 13, 0, 22, 0, 56,
        0, 10, 6, 6, 110, 111, 110, 97,109, 101};
    int gdspostamble[8]     = {0, 4, 7, 0, 0, 4, 4, 0};
    unsigned char gdsPre[102];
    for (int k = 0; k < 102; k++)
        gdsPre[k] = (unsigned char) gdspreamble[k];
    for (int k = 0; k < 8; k++)
        gdsPost[k] = (unsigned char) gdspostamble[k];
    fwrite(gdsPre, sizeof(char), 96, outputFile);
    fwrite(filename, sizeof(char), 6, outputFile);
}

void readcor(FILE * outputFile,int p,int q){
    FILE *tempFile= NULL;
    char filename[20]="test";
    char *buffer;
    long lSize;
    char *cCoords;
    size_t result;
    strcat(strcat(strcat(strcat(filename,num2str(p)),"&"),num2str(q)),".cor");
    if ((tempFile = fopen(filename, "rb")) == NULL)
        printf("Cannot open file.\n");
    fseek (tempFile, 0 , SEEK_END);
    lSize = ftell (tempFile);
    rewind (tempFile);
    buffer = (char*) malloc (lSize);
    result=fread(buffer, sizeof(char), lSize, tempFile);
    cCoords=new char [lSize];
//     encodePoly32(buffer, cCoords,lSize/4);
    fwrite(buffer, sizeof(char), lSize, outputFile);
//      std::cout << std::fixed <<'n'<< '\t'<<(int)'n'<< '\n';
    fclose (tempFile);
    free (buffer);
    delete cCoords;
}


// Gateway function for MATLAB interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int Mx      = mxGetScalar(prhs[0]);
    int My      = mxGetScalar(prhs[1]);
    int clen    =  1+ mxGetN(prhs[2]);
    char * filename;
    filename= (char*) malloc(clen);
    mxGetString(prhs[2], filename, (mwSize)clen);
    FILE * outputFile = NULL;
    unsigned char gdsPost[8];
    if ((outputFile = fopen(filename, "wb")) == NULL)
                printf("Cannot open file.\n");
    initGDS(outputFile, gdsPost,filename,clen);
    for (int p=1;p<=Mx;p++){
        for (int q=1;q<=My;q++){
            readcor(outputFile,p,q);
        }
    }
    renderGDS(outputFile, gdsPost);   
}
