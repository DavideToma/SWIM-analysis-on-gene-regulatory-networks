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

function [] = connectedComponents_plot(soglia,thr,ratio,dirnameFigureSwitch)

figure('Color', 'w')

plot(thr,ratio,'b-','LineWidth',2)
hold on
plot([soglia soglia],[min(ratio) max(ratio)],'r--','LineWidth',2)
ylim([min(ratio) max(ratio)+0.1])
xlabel('Correlation threshold',...
    'fontsize',12,'fontname','Times New Roman', 'FontWeight','bold')
ylabel('Fraction of nodes of the largest component','fontsize',12,'fontname',...
    'Times New Roman', 'FontWeight','bold')

saveas(gcf,strcat(dirnameFigureSwitch,'connectedComponent.fig'))
saveas(gcf,strcat(dirnameFigureSwitch,'connectedComponent.png'))

end
