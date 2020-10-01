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

function []=trimodal(dirnameFigureSwitch,APCC)

hubs=cell2mat(APCC(:,1));
ind=strmatch('no hub',APCC(:,3),'exact');
hubs(ind)=[];

figure('Color','w')

figure(1)
[f0,x0] = ksdensity(hubs); 
plot(x0,f0,'-k');
xlabel('Average Pearson Correlation Coefficient (APCC)',...
    'fontsize',12,'fontname','Times New Roman', 'FontWeight','bold')
ylabel('Probability density','fontsize',12,'fontname',...
    'Times New Roman', 'FontWeight','bold')

saveas(gcf,strcat(dirnameFigureSwitch,'APCC-distribution.fig'))
saveas(gcf,strcat(dirnameFigureSwitch,'APCC-distribution.png'))

end







