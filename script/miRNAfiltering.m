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

function [CpzmiRNA,NpzmiRNA,NmiRNA,CmiRNA,miRNAid]=...
    miRNAfiltering(CpzmiRNA,NpzmiRNA,NmiRNA,CmiRNA,miRNAid,meanCmir,meanNmir,...
    fc_soglia,FDR_soglia,ok_log,dirnameTxtFilter)

pvalues = mattest(NmiRNA, CmiRNA);
fdr = mafdr(pvalues, 'BHFDR', true);

if(ok_log)
    mavolcanoplot(CmiRNA,NmiRNA,fdr,'Labels',miRNAid,'PCUTOFF', FDR_soglia,...
        'FOLDCHANGE', fc_soglia);
else
    mavolcanoplot(CmiRNA,NmiRNA,fdr,'Labels',miRNAid,'PCUTOFF', FDR_soglia,...
        'FOLDCHANGE', fc_soglia,'LogTrans','true');
end

ind = find(fdr > FDR_soglia);

miRNAid(ind)=[];
NmiRNA(ind,:)=[];
CmiRNA(ind,:)=[];
pvalues(ind)=[];
fdr(ind)=[];
meanCmir(ind)=[];
meanNmir(ind)=[];

if(ok_log)
    logFC = mean(CmiRNA,2) - mean(NmiRNA,2);
    
    figure
    [freq,xcenter]=hist(logFC,100);
    h=bar(xcenter,freq,'hist');
    
    x1 = abs(xcenter-log2(fc_soglia)) < diff(xcenter(1:2))/2; % Find the index of the bin containing the threshold
    x2 = abs(xcenter+log2(fc_soglia)) < diff(xcenter(1:2))/2; % Find the index of the bin containing the threshold
    
    index=(xcenter>=xcenter(x1==1) | xcenter<=xcenter(x2==1));
    colors = [(index(:)+0.5.*(~index(:))) 0.5.*(~index(:)) 0.5.*(~index(:))];
    set(h,'FaceVertexCData',colors);
    xlabel('Log2 (ratio)',...
        'fontsize',12,'fontname','Times New Roman', 'FontWeight','bold')
    ylabel('Frequency','fontsize',12,'fontname',...
        'Times New Roman', 'FontWeight','bold')
    title('microRNA','fontsize',12,'fontname',...
        'Times New Roman', 'FontWeight','bold')
    xlim([min(logFC) max(logFC)])
    
    ind = ( abs(logFC) <= log2(fc_soglia) );
    
    miRNAid(ind)=[];
    NmiRNA(ind,:)=[];
    CmiRNA(ind,:)=[];
    pvalues(ind)=[];
    fdr(ind)=[];
    logFC(ind)=[];
    meanCmir(ind)=[];
    meanNmir(ind)=[];

    data=[NmiRNA CmiRNA];
    pz=[NpzmiRNA; CpzmiRNA];
    N=length(miRNAid);
    
    filename = strcat(dirnameTxtFilter,'stat-log-dataFiltered-miRNA.txt');
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', 'gene', 'pvalues', 'fdr', ...
        'mean cancer', 'mean normal', 'log-ratio');
    for i = 1 : N
        fprintf(fid, '%s\t%f\t%f\t%f\t%f\t%f\n', miRNAid{i}, pvalues(i), fdr(i), ...
            meanCmir(i), meanNmir(i), logFC(i));
    end
    fclose(fid);
    
    filename = strcat(dirnameTxtFilter,'log-dataFiltered-miRNA.txt');
    fid = fopen(filename, 'w');
    fprintf(fid, 'gene');
    for c = 1 : numel(pz)
        fprintf(fid, '\t%s', char(pz(c)));
    end
    fprintf(fid, '\n');
    for i = 1 : N
        fprintf(fid, '%s', miRNAid{i});
        for c = 1 : numel(pz)
            fprintf(fid, '\t%f', data(i,c));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
else
    disp('------------------------------------------------------')
    disp('ATTENZIONE: non hai fatto il log2 dei dati di microRNA')
    disp('------------------------------------------------------')
    
    FC = ( meanCmir + 1) ./ ( meanNmir + 1);
    
    figure
    [freq,xcenter]=hist(FC,100);
    h=bar(xcenter,freq,'hist');
    
    x1 = abs(xcenter-fc_soglia) < diff(xcenter(1:2))/2; % Find the index of the bin containing the threshold
    x2 = abs(xcenter-(1/fc_soglia)) < diff(xcenter(1:2))/2; % Find the index of the bin containing the threshold
    
    if(any(x1) && any(x2))
        index = (xcenter>=xcenter(x1==1) | xcenter<=xcenter(x2==1));
        colors = [(index(:)+0.5.*(~index(:))) 0.5.*(~index(:)) 0.5.*(~index(:))];
    else
        colors = ones(numel(xcenter),1) * [0.5 0.5 0.5];
    end
    set(h,'FaceVertexCData',colors);
    xlabel('Ratio',...
        'fontsize',12,'fontname','Times New Roman', 'FontWeight','bold')
    ylabel('Frequency','fontsize',12,'fontname',...
        'Times New Roman', 'FontWeight','bold')
    title('microRNA','fontsize',12,'fontname',...
        'Times New Roman', 'FontWeight','bold')
    xlim([0 max(FC)])
    
    ind = find( (FC <= fc_soglia) & (FC >= (1/fc_soglia)) );
    
    miRNAid(ind)=[];
    NmiRNA(ind,:)=[];
    CmiRNA(ind,:)=[];
    pvalues(ind)=[];
    fdr(ind)=[];
    FC(ind)=[];
    meanCmir(ind)=[];
    meanNmir(ind)=[];
    
    data=[NmiRNA CmiRNA];
    pz=[NpzmiRNA; CpzmiRNA];
    N=length(miRNAid);
    
    filename = strcat(dirnameTxtFilter,'stat-dataFiltered-miRNA.txt');
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', 'gene', 'pvalues', 'fdr', ...
        'mean cancer', 'mean normal', 'ratio');
    for i = 1 : N
        fprintf(fid, '%s\t%f\t%f\t%f\t%f\t%f\n', miRNAid{i}, pvalues(i), fdr(i), ...
            meanCmir(i), meanNmir(i), FC(i));
    end
    fclose(fid);
    
    filename = strcat(dirnameTxtFilter,'dataFiltered-miRNA.txt');
    fid = fopen(filename, 'w');
    fprintf(fid, 'gene');
    for c = 1 : numel(pz)
        fprintf(fid, '\t%s', char(pz(c)));
    end
    fprintf(fid, '\n');
    for i = 1 : N
        fprintf(fid, '%s', miRNAid{i});
        for c = 1 : numel(pz)
            fprintf(fid, '\t%f', data(i,c));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end

end