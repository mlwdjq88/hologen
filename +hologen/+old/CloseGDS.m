function CloseGDS(outputFile)
gdsPost = [0, 4, 7, 0, 0, 4, 4, 0];
fwrite(outputFile,gdsPost , 'uint8' );
fclose(outputFile);