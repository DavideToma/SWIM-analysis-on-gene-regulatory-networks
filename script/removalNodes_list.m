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

function [nodeTOremove_list] = removalNodes_list(N,nodeName0,nodeName,...
    ind_list,degree_list)

if(~isempty(degree_list))
    [val,ind]=sort(degree_list,'descend');
end

if(~isempty(ind_list))
    list=nodeName(ind_list(ind));
elseif(~isempty(degree_list))
    list=nodeName(ind);
else
    list=nodeName;
end

M=numel(list);

if(N<=M)
    lim=N;
else
    lim=M;
end

nodeTOremove_list=zeros(lim,1);
for i=1:lim
    nodeTOremove_list(i)=strmatch(list(i),nodeName0,'exact');
end

end

