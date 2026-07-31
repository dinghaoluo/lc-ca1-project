% SaveRes(FileName, Res)
%
% Saves a .res file from an array -

function SaveRes(FileName, Res)

outputfile = fopen(FileName,'w');
fprintf(outputfile,'%d\n', Res(:));
fclose(outputfile);
