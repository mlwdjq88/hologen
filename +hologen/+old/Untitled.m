exportPolygon(trapCoords, polyPre, polyPost, polyForm, outputFile, File_format, clockSpeed);
initGDS(outputFile, gdsPost, polyPre, polyPost, polyForm, layerNumber);


void initGDS(FILE * outputFile, unsigned char * gdsPost, unsigned char* polyPre, unsigned char * polyPost, unsigned char * polyForm, int layerNumber){
    int gdspreamble[102] = { 0, 6, 0, 2, 0, 7, 0, 28, 1, 2, 230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
        230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 10, 2, 6, 110, 111, 110, 97,
        109, 101, 0, 20, 3, 5, 61, 104, 219, 139, 172, 113, 12, 180, 56, 109,
        243, 127, 103, 94, 246, 236, 0, 28, 5, 2, 0, 114, 0, 4, 0, 17, 0, 13,
        0, 22, 0, 56, 0, 114, 0, 4, 0, 17, 0, 13, 0, 22, 0, 56,
        0, 10, 6, 6, 110, 111, 110, 97, 109, 101};
    
    int gdspostamble[8]     = {0, 4, 7, 0, 0, 4, 4, 0};
    
    // {0, 4, 8, 0, 0, 6, 13, 2, 0, [layer = 1], 0, 6, 14, 2, 0, 0};
    int polypreamble[16]    = {0, 4, 8, 0, 0, 6, 13, 2, 0, layerNumber, 0, 6, 14, 2, 0, 0};
    int polypostamble[4]    = {0, 4, 17, 0};
    int polyBlockFormat[4]  = {0, 44, 16, 3};
    
    unsigned char gdsPre[102];
    for (int k = 0; k < 102; k++)
        gdsPre[k] = (unsigned char) gdspreamble[k];
    for (int k = 0; k < 8; k++)
        gdsPost[k] = (unsigned char) gdspostamble[k];
    for (int k = 0; k < 16; k++)
        polyPre[k] = (unsigned char) polypreamble[k];
    for (int k = 0; k < 4; k++)
        polyPost[k] = (unsigned char) polypostamble[k];
    for (int k = 0; k < 4; k++)
        polyForm[k] = (unsigned char) polyBlockFormat[k];

    fwrite(gdsPre, sizeof(char), 102, outputFile);
}
void exportPolygon(long * coords, unsigned char * polyPre, unsigned char * polyPost,
                   unsigned char * polyForm, FILE * outputFile, int File_format, long clockSpeed){
    switch (File_format){
        case 1:
        case 2:
            char cCoords[40];
            fwrite(polyPre, sizeof(char), 16, outputFile);
            fwrite(polyForm, sizeof(char), 4, outputFile);
            encodePoly32(coords, cCoords);
            fwrite(cCoords, sizeof(char), 40, outputFile);
            fwrite(polyPost, sizeof(char), 4, outputFile);
            break;
        case 3:
            fractureAndWriteWRVPoly(coords, outputFile, clockSpeed);
            break;
        }
}

void encodePoly32(long * coords, char * cCoords){
    char cPart[4];
    for (int k = 0; k < 10; k++){
        encode32(coords[k], cPart);
        for (int m = 0; m < 4; m++){
            cCoords[k*4 + m] = cPart[m];
        }
    }
}

void encode32(long aCoord, char * cPart){
    cPart[0] = (aCoord >> 24) & 255;
    cPart[1] = (aCoord >> 16) & 255;
    cPart[2] = (aCoord >> 8) & 255;
    cPart[3] = (aCoord) & 255;
}
