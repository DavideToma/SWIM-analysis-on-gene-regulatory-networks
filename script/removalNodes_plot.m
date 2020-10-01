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

function [] = removalNodes_plot(ind_date,ind_party,fH,lambdaH,...
    fFC,lambdaFC,fD,lambdaD,fP,lambdaP,fR,lambdaR,fNS,lambdaNS,...
    fS,lambdaS,dirnameFigureSwitch)

figure('Color','w')

subplot(1,2,1)
plot(fH,lambdaH,'m-','LineWidth',2)
ylim([(min(lambdaR)-0.01) (max(lambdaH)+0.05)])
xlabel('Fraction of removed nodes',...
    'fontsize',12,'fontname','Times New Roman', 'FontWeight','bold')
ylabel('Average shortest path','fontsize',12,'fontname',...
    'Times New Roman', 'FontWeight','bold')
hold on
plot(fFC,lambdaFC,'r-','LineWidth',2)

if(~isempty(ind_date))
    plot(fD,lambdaD,'k-','LineWidth',2)
end

if(~isempty(ind_party))
    plot(fP,lambdaP,'b','LineWidth',2)
end
plot(fR,lambdaR,'g-','LineWidth',2)

if(~isempty(ind_date) && ~isempty(ind_party))
    legend('hubs','fight-club', 'date','party','random')
elseif(isempty(ind_date)  && ~isempty(ind_party))
    legend('hubs','fight-club','party','random')
elseif(~isempty(ind_date)  && isempty(ind_party))
    legend('hubs','fight-club', 'date','random')
end

subplot(1,2,2)
plot(fNS,lambdaNS,'b-','LineWidth',2)
ylim([(min(lambdaR)-0.01) (max(lambdaH)+0.05)])
xlabel('Fraction of removed nodes',...
    'fontsize',12,'fontname','Times New Roman', 'FontWeight','bold')
ylabel('Average shortest path','fontsize',12,'fontname',...
    'Times New Roman', 'FontWeight','bold')
hold on
plot(fS,lambdaS,'r-','LineWidth',2)
plot(fR,lambdaR,'g-','LineWidth',2)
legend('non-switch hubs','switch','random')

saveas(gcf,strcat(dirnameFigureSwitch,'removalNodes.fig'))
saveas(gcf,strcat(dirnameFigureSwitch,'removalNodes.png'))

end

