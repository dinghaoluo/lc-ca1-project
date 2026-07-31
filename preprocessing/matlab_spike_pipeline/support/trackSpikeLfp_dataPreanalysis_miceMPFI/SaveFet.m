% SaveFet(FileName, Fet)
%
% Saves a feature file

function SaveFet(FileName, Fet);

nFeatures = size(Fet, 2);
formatstring = '%d';
for ii=2:nFeatures
  formatstring = [formatstring,'\t%d'];
end
formatstring = [formatstring,'\n'];

outputfile = fopen(FileName,'w');
fprintf(outputfile, '%d\n', nFeatures);
fprintf(outputfile,formatstring,round(Fet'));
fclose(outputfile);
