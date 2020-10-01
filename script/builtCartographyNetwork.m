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

function [adjMatrix,nodeName] = builtCartographyNetwork(rho,diag_Nr,diag_Nc,net,...
    nodeName,adjMatrix,nodeNameG,dirnameMatSwitch,dirnameTxtSwitch)

[name,ind]=setdiff(nodeName,nodeNameG');
nodeName(ind)=[];
adjMatrix(ind,:)=[];
adjMatrix(:,ind)=[];

count=1;
while(count)
    
    [name1,ind_name1,ind_net1]=intersect(name,net(:,1));
    [name2,ind_name2,ind_net2]=intersect(name,net(:,2));
    
    if(~isempty(ind_net1))
        diag_Nr(ind_net1)=[];
        diag_Nc(ind_net1)=[];
    end
    
    if(~isempty(ind_net2))
        diag_Nr(ind_net2)=[];
        diag_Nc(ind_net2)=[];
    end
    
    net(unique([ind_net1;ind_net2]),:)=[];
    
    if ( isempty(ind_net1) && isempty(ind_net2))
        count=0;
    end
end

save(strcat(dirnameMatSwitch,'CartographyNetwork'),'nodeName','net')

fid=fopen(strcat(dirnameTxtSwitch,'CartographyNetwork.txt'),'w');
fprintf(fid,'%s\t%s\t%s\n','Source','Target', 'Correlation');
for i=1:size(net,1)
    fprintf(fid,'%s\t%s\t%f\n',char(net(i,1)),char(net(i,2)),...
        rho(diag_Nr(i),diag_Nc(i)));
end
fclose(fid);

end

