% Copyright 2016 Paola Paci
%
% This file is part of SWIM.
%
% SWIM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SWIM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SWIM.  If not, see <http://www.gnu.org/licenses/>.
%

function [CpzmiRNA,NpzmiRNA,miRNAid,CmiRNA,NmiRNA,meanCmir,meanNmir]=...
    miRNAextractData(dirnameMatFilter,dirnameTxtFilter,fnamePzC,fnamePzN,fnameData)

tmp=importdata(fnameData);
miRNAid=tmp.textdata(2:end,1);

CpzmiRNA=importdata(fnamePzC);
NpzmiRNA=importdata(fnamePzN);

CmiRNA=selectPz(CpzmiRNA,tmp);
NmiRNA=selectPz(NpzmiRNA,tmp);

meanCmir=mean(CmiRNA,2);
meanNmir=mean(NmiRNA,2);

pvalues = mattest(NmiRNA, CmiRNA);
fdr = mafdr(pvalues, 'BHFDR', true);

save(strcat(dirnameMatFilter,'/dataORIG-miRNA'),'CpzmiRNA',...
    'NpzmiRNA','miRNAid','CmiRNA','NmiRNA','meanCmir','meanNmir')

N=numel(miRNAid);
filename = strcat(dirnameTxtFilter,'/stat-dataORIG-miRNA.txt');
fid = fopen(filename, 'w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', 'gene', 'pvalues', 'fdr', ...
    'mean cancer', 'mean normal');
for i = 1 : N
    fprintf(fid, '%s\t%f\t%f\t%f\t%f\n', miRNAid{i}, pvalues(i), fdr(i), ...
        meanCmir(i), meanNmir(i));
end
fclose(fid);

end


