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

function idx=kmean_II(MaxIter,num_rep,num_cluster,distMatrix,nodeName,...
    dirnameMatSwitch,dirnameTxtSwitch)

opts=statset('Display','final','MaxIter',MaxIter);
fprintf('Number of clusters = %d\n', num_cluster);
[idx,c,sumd]=kmeans(distMatrix,num_cluster,'replicates',num_rep,'options',opts);

M=numel(nodeName);

save(strcat(dirnameMatSwitch,'idx'),'idx')

filename=strcat(dirnameTxtSwitch,'idx.txt');
fid=fopen(filename,'w');
for i=1:M
   fprintf(fid,'%s\t%d\n',char(nodeName(i)),idx(i));
end
fclose(fid);

end

