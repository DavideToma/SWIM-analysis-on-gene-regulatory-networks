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

function [] = connectedComponents(rho,soglia,minimum,step,maximum,...
    dirnameMatSwitch,dirnameFigureSwitch)

thr=minimum:step:maximum;

ratio=zeros(numel(thr),1);
count=0;
for h = minimum : step : maximum
    
    [Nr,Nc]=find(abs(rho)>=h);
    
    nodi=unique(Nr); 
    
    M=numel(nodi);
    adjMatrix=zeros(M);
    for i=1:M
        found=find(nodi(i)==Nr);
        c=numel(found);
        for k=1:c
            j=find(Nc(found(k))==nodi);
            adjMatrix(i,j)=1;
        end
    end
    
    sparse_adj=sparse(adjMatrix);
    
    [S,C]=graphconncomp(sparse_adj,'Directed',false);
    tmp=tabulate(C);
    numel_S=sort(tmp(:,2),'descend');
    numel_Smax=numel_S(1);
   
    count=count+1;
    ratio(count)=numel_Smax/M;
       
end

save(strcat(dirnameMatSwitch,'connectedComponent'), 'thr', 'ratio')

connectedComponents_plot(soglia,thr,ratio,dirnameFigureSwitch)

end
