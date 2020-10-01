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

function data=selectPz(pz,tmp)

list = tmp.textdata(1,2:end);

N = numel(pz);
M = size(tmp.data,1);

data = zeros(M,N);
for i = 1 : N
    found = strmatch(pz{i},list,'exact');
    data(:,i) = tmp.data(:,found);
end

end

