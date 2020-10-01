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

function [f,size,lambda] = connectivityChange(N,size0,nodeTOremove_list,...
    adj0,nodeName0)

M=numel(nodeTOremove_list);

if(N<=M)
    lim=N;
else
    lim=M;
end
    
f=[];
size=[];
lambda=[];
nodeTOremove=[];
for i=1:lim
    
    adj_tmp=adj0;
    nodeName_tmp=nodeName0;
    
    nodeTOremove=[nodeTOremove nodeTOremove_list(i)];
    
    adj_tmp(nodeTOremove,:)=[];
    adj_tmp(:,nodeTOremove)=[];
    nodeName_tmp(nodeTOremove)=[];
    
    [size_tmp adj nodeName degree_tmp lambda_tmp]=connectivity(adj_tmp,...
        nodeName_tmp);
     
    size(i)=size_tmp;
    lambda(i)=lambda_tmp;
    f(i)=numel(nodeTOremove)/size0;
end

end

