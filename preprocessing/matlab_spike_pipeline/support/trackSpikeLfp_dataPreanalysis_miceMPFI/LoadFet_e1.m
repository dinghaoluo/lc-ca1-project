% Fet = LoadFet(FileName)
%
% A simple matlab function to load a .fet file

function Fet = LoadFet_e1(FileName);

Fp = fopen(FileName, 'r');

if Fp==-1
    error(['Could not open file ' FileName]);
end

nFeatures = fscanf(Fp, '%d', 1);
Fet = fscanf(Fp, '%f', [nFeatures, inf])';

fclose(Fp);

return;