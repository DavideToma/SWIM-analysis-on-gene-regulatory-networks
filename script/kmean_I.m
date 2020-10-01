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

function [distMatrix,adjMatrix]=kmean_I(Nr,Nc,nodi,rho,MaxIter,num_rep,...
    N,dirnameMatSwitch)

M=numel(nodi);
distMatrix=zeros(M);
adjMatrix=zeros(M);
for i=1:M
    found=find(nodi(i)==Nr);
    c=numel(found);
    for k=1:c
        j=find(Nc(found(k))==nodi);
        distMatrix(i,j)=1-rho(Nr(found(k)),Nc(found(k)));
        adjMatrix(i,j)=1;
    end
end

opts=statset('Display','final','MaxIter',MaxIter);
SSE=[];
for i=1:N
    fprintf('Number of clusters = %d\n', i);
    [idx,c,sumd]=kmeans(distMatrix,i,'replicates',num_rep,'options',opts);
    SSE(i)=sum(sumd);
end

save(strcat(dirnameMatSwitch,'SSE'),'SSE')
save(strcat(dirnameMatSwitch,'adjMatrix'),'adjMatrix')

figure('Color','w')
plot(SSE,'o-')
xlabel('Number of clusters','fontsize',12,'fontname',...
    'Times New Roman', 'FontWeight','bold')
ylabel('Sum of the Squared Error (SSE)','fontsize',12,'fontname',...
    'Times New Roman', 'FontWeight','bold')
title('Scree plot','fontsize',12,'fontname',...
    'Times New Roman', 'FontWeight','bold')

end
