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

function [CpzmRNA,NpzmRNA,mRNAid,CmRNA,NmRNA,meanC,meanN]=...
    mRNAextractData(dirnameMatFilter,dirnameTxtFilter,fnamePzC,fnamePzN,fnameData)

tmp=importdata(fnameData);
mRNAid=tmp.textdata(2:end,1);

CpzmRNA=importdata(fnamePzC);
NpzmRNA=importdata(fnamePzN);

CmRNA=selectPz(CpzmRNA,tmp);
NmRNA=selectPz(NpzmRNA,tmp);

meanC=mean(CmRNA,2);
meanN=mean(NmRNA,2);

pvalues = mattest(NmRNA, CmRNA);
fdr = mafdr(pvalues, 'BHFDR', true);

save(strcat(dirnameMatFilter,'/dataORIG'),'CpzmRNA','NpzmRNA',...
    'CmRNA','NmRNA','mRNAid','meanC','meanN')

N=numel(mRNAid);
filename = strcat(dirnameTxtFilter,'/stat-dataORIG.txt');
fid = fopen(filename, 'w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', 'gene', 'pvalues', 'fdr', ...
    'mean cancer', 'mean normal');
for i = 1 : N
    fprintf(fid, '%s\t%f\t%f\t%f\t%f\n', mRNAid{i}, pvalues(i), fdr(i), ...
        meanC(i), meanN(i));
end
fclose(fid);

end


