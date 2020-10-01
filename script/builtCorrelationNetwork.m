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

function [soglia,Nr,Nc,diag_Nr,diag_Nc,nodi,nodeName,net]=builtCorrelationNetwork(rho,percentile,...
    RNAname,dirnameTxtSwitch,dirnameMatSwitch)

soglia=prctile(rho(:),percentile);

[freq,xcenter]=hist(rho(:),100);
h=bar(xcenter,freq,'hist');

x1 = abs(xcenter-soglia) < diff(xcenter(1:2))/2; % Find the index of the bin containing the threshold
x2 = abs(xcenter+soglia) < diff(xcenter(1:2))/2; % Find the index of the bin containing the threshold

index=(xcenter>=xcenter(x1==1) | xcenter<=xcenter(x2==1));
colors = [(index(:)+0.5.*(~index(:))) 0.5.*(~index(:)) 0.5.*(~index(:))]; 
set(h,'FaceVertexCData',colors);
xlabel('Pearson correlation coefficient',...
    'fontsize',12,'fontname','Times New Roman', 'FontWeight','bold')
ylabel('Frequency','fontsize',12,'fontname',...
    'Times New Roman', 'FontWeight','bold')
xlim([-1 1])

% hold on
% y1 = freq(x1);
% y2 = freq(x2);
% plot([soglia soglia], [0 y1],'k--','LineWidth',1.5)
% plot([-soglia -soglia], [0 y2],'k--','LineWidth',1.5)
% hold off

[Nr, Nc]=find(abs(rho)>=soglia);
nodi=unique(Nr);                 % the same as unique([Nr;Nc])
nodeName=RNAname(nodi);

diag=tril(rho);
[diag_Nr, diag_Nc]=find(abs(diag)>=soglia);

net=[RNAname(diag_Nr), RNAname(diag_Nc)];

save(strcat(dirnameMatSwitch,'CorrelationNetwork'), 'nodeName', 'net',...
    'diag_Nr','diag_Nc')

fid=fopen(strcat(dirnameTxtSwitch,'CorrelationNetwork.txt'),'w');
fprintf(fid,'%s\t%s\t%s\n','Source','Target', 'Correlation');
for i=1:size(net,1)
    fprintf(fid,'%s\t%s\t%f\n',char(net(i,1)),char(net(i,2)),...
        rho(diag_Nr(i),diag_Nc(i)));
end
fclose(fid);

end