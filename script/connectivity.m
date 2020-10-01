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

function [max_com_size adj_max_conn_com nodeName_max_conn_com degree_max_conn_com lambda]=connectivity(adjMatrix, nodeName)

G = sparse(adjMatrix);
[S C]=graphconncomp(G,'Directed',false);

node_S=cell(1,S);
size_S=[];
for i=1:S
    ind=find(C==i);
    node_S{1,i}=ind';
    size_S(i)=numel(ind);
end
[max_com_size ind_max_com_size]=max(size_S);
nodeRow_max_conn_com=node_S{1,ind_max_com_size};
nodeName_max_conn_com=nodeName(nodeRow_max_conn_com);

adj_max_conn_com=adjMatrix(nodeRow_max_conn_com,nodeRow_max_conn_com);
G = sparse(adj_max_conn_com);
D = graphallshortestpaths(G);
% http://code.google.com/p/visualconnectome/source/browse/trunk/Plugins/BCT/charpath.m?r=2
% average shortest path
lambda = sum(sum(D(D~=Inf)))/length(nonzeros(D~=Inf));

degree_max_conn_com=zeros(max_com_size,1);
for i=1:max_com_size
    degree_max_conn_com(i)=sum(nonzeros(adj_max_conn_com(i,:)>0));
end

end
