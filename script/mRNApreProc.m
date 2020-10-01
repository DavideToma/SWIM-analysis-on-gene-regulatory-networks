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

function [CpzmRNA,NpzmRNA,NmRNA,CmRNA,mRNAid,meanC,meanN]=...
    mRNApreProc(CpzmRNA,NpzmRNA,NmRNA,CmRNA,mRNAid,meanC,meanN,...
    ok_log,filtro_zeri,filtro_iqr)

if(ok_log)
    NmRNA = log2(NmRNA + 1);
    CmRNA = log2(CmRNA + 1);
    data = [NmRNA CmRNA];
    iqr_list = iqr(data');
else
    data = [NmRNA CmRNA];
    iqr_list = log2(iqr(data') + 1);
end

soglia_iqr = prctile(iqr_list,filtro_iqr);

figure(1)

[freq,xcenter]=hist(iqr_list,100);
h=bar(xcenter,freq,'hist');
if(ok_log)
    xlabel('Interquartile range (IQR)',...
        'fontsize',12,'fontname','Times New Roman', 'FontWeight','bold')
else
    xlabel('Log2 of the interquartile range (IQR)',...
        'fontsize',12,'fontname','Times New Roman', 'FontWeight','bold')
end
ylabel('Frequency','fontsize',12,'fontname',...
    'Times New Roman', 'FontWeight','bold')
title('RNA','fontsize',12,'fontname',...
    'Times New Roman', 'FontWeight','bold')
xlim([0 max(iqr_list)])
    
x = abs(xcenter-soglia_iqr) <= diff(xcenter(1:2))/2; % Find the index of the bin containing the threshold

if(any(x))
    index=(xcenter<=xcenter(x==1));
    colors = [(~index(:)+0.5.*(index(:))) 0.5.*(index(:)) 0.5.*(index(:))];
else
    colors = ones(numel(xcenter),1)*[1 0 0];
end

set(h,'FaceVertexCData',colors);

M = numel(mRNAid);
N = numel(CpzmRNA) + numel(NpzmRNA);
prc_zeri = zeros(M,1);
for i = 1 : M
    found = find(data(i,:)==0);
    prc_zeri(i) = numel(found) * 100 / N;
end

figure(2)
scatter(iqr_list,100-prc_zeri,'MarkerEdgeColor',[0.5 0.5 0.5],...
    'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1.)
title('RNA','fontsize',12,'fontname','Times New Roman',...
    'FontWeight','bold')
if(ok_log)
    xlabel('Interquartile range (IQR)','fontsize',12,...
        'fontname','Times New Roman', 'FontWeight','bold')
else
    xlabel('Log2 of the interquartile range (IQR)','fontsize',12,...
        'fontname','Times New Roman', 'FontWeight','bold')
end
ylabel('100 - % of zeros','fontsize',12,'fontname','Times New Roman',...
    'FontWeight','bold')
hold on
plot([0 max(iqr_list)], [100-filtro_zeri 100-filtro_zeri],'k--','LineWidth',1.5)
plot([soglia_iqr soglia_iqr], [0 100],'k--','LineWidth',1.5)
xlim([0 max(iqr_list)])

ind = find(iqr_list <= soglia_iqr);

if(~isempty(ind))
    mRNAid(ind)=[];
    NmRNA(ind,:)=[];
    CmRNA(ind,:)=[];
    prc_zeri(ind)=[];
    iqr_list(ind)=[];
    meanC(ind)=[];
    meanN(ind)=[];
end

clear ind

ind = find(prc_zeri > filtro_zeri);
        
if(~isempty(ind))
    mRNAid(ind)=[];
    NmRNA(ind,:)=[];
    CmRNA(ind,:)=[];
    prc_zeri(ind)=[];
    iqr_list(ind)=[];
    meanC(ind)=[];
    meanN(ind)=[];
end

scatter(iqr_list,100-prc_zeri,'MarkerEdgeColor',[1 0 0],...
    'MarkerFaceColor',[1 0 0])
hold off

end
          


