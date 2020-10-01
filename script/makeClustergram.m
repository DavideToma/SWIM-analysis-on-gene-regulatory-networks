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

function [] = makeClustergram(data,RNAname,tissue,CpzmRNA,NpzmRNA,...
    dirname,name,title)

colour_red = {};
for i = 1 : numel(CpzmRNA)
    colour_red{i} = 'r';
end
colour_black = {};
for i = 1 : numel(NpzmRNA)
    colour_black{i} = 'k';
end
color_tissue = [colour_red'; colour_black'];

C = exp_colormap('blue-yellow',64);

CGobj = clustergram(data,'RowLabels',RNAname,'ColumnLabels',tissue,...
    'Standardize','row','Cluster','all','RowPdist','correlation',...
    'ColumnPdist','correlation','Linkage','complete',...
    'DisplayRange',3,'Symmetric','true', 'Colormap',C,...
    'ColumnLabelsColor',struct('Labels', tissue, 'Colors', color_tissue));
addTitle(CGobj, title)
fig = plot(CGobj);
saveas(fig,strcat(dirname,strcat(name,'.fig')))
saveas(fig,strcat(dirname,strcat(name,'.png')))
delete(CGobj)
close all


end

