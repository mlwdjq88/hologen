function exportPolygon(outputFile,coords, polyPre, polyPost,polyForm )
            fwrite(outputFile,polyPre, 'uint8' );
            fwrite(outputFile,polyForm,'uint8' );
%             for i=1:10
%                 str=dec2hex(coords(i),8);
%                 fwrite(outputFile,hex2dec(str(1:2)), 'uint8');
%                 fwrite(outputFile,hex2dec(str(3:4)), 'uint8');
%                 fwrite(outputFile,hex2dec(str(5:6)), 'uint8');
%                 fwrite(outputFile,hex2dec(str(7:8)), 'uint8');
% %             
%             end
             fwrite(outputFile,coords', 'int32','b');
            fwrite(outputFile,polyPost,'uint8');

