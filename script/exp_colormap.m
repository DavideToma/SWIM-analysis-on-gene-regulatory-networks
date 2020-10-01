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

function C=exp_colormap(colors, N)
% a simple function that returns a colormap, C, for visualizing gene
% expression.  C is just a N x 3 matrix [R G B] describing the range of color values.
%
% example usage:
%    C = exp_colormap('blue-yellow',64);
%    colormap(C);
% 
% called without any arguments, it returns a [3 x 64] green-red colormap.
%
% options for 'colors' are:
%
%  'blue-yellow'
%  'green-red'
%  'yellow'
%
% 'N' represents the number of degrees of color to use.  the default is 64.
%
% generally speaking, the two-colored maps are appropriate for visualizing
% expression data which is normalized to a reference condition, e.g. to show
% whether expression is lower (blue or green) or higher (yellow or red) than
% the reference.
%
% the single-color yellow map ('yellow') is appropriate for displaying
% levels of gene expression which are not compared (necessarily) to a single
% reference, and this is similar to the colormap used in the D-chip
% expression analysis software.
% 

% the colormaps returned range monotonically.
%  
  
if 1 ~= exist('colors') 
  colors = 'green-red';
end

if 1 ~= exist('N') 
  N = 64;
end

X = [0.5: -1/(N-1):-0.5];
X = abs(X).*2;

switch colors
case {'green-red'}
  R = [X(:, 1:N/2) zeros(1,N/2)];
  G = [zeros(1, N/2) X(:,(N/2 + 1):N)];
  B = zeros(1,N);

case {'blue-yellow'} 
  R = [zeros(1,N/2) X(:,(N/2 + 1):N)];
  B = [X(:,1:N/2) zeros(1,N/2)];
  G = [zeros(1,N/2) X(:,(N/2 + 1):N)];

case {'yellow'} 
  X = [0:1/(N - 1):1];
  R = X;
  B = zeros(1,N);
  G = X;

otherwise
 error([colors ' is not a known option for coloring']);
end

C = [R' G' B'];
