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

function [CpzmRNA,NpzmRNA,NmRNA,CmRNA,mRNAid]=...
    mRNAfiltering(CpzmRNA,NpzmRNA,NmRNA,CmRNA,mRNAid,meanC,meanN,...
    fc_soglia,FDR_soglia,ok_log,dirnameTxtFilter)

pvalues = mattest(NmRNA, CmRNA);
fdr = mafdr(pvalues, 'BHFDR', true);

if(ok_log)
    mavolcanoplot(CmRNA,NmRNA,fdr,'Labels',mRNAid,'PCUTOFF', FDR_soglia,...
        'FOLDCHANGE', fc_soglia);
else
    mavolcanoplot(CmRNA,NmRNA,fdr,'Labels',mRNAid,'PCUTOFF', FDR_soglia,...
        'FOLDCHANGE', fc_soglia,'LogTrans','true');
end

ind = find(fdr > FDR_soglia);

mRNAid(ind)=[];
NmRNA(ind,:)=[];
CmRNA(ind,:)=[];
pvalues(ind)=[];
fdr(ind)=[];
meanC(ind)=[];
meanN(ind)=[];

if(ok_log)
    logFC = mean(CmRNA,2) - mean(NmRNA,2);
    
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
    title('RNA','fontsize',12,'fontname',...
        'Times New Roman', 'FontWeight','bold')
    xlim([min(logFC) max(logFC)])

    ind = ( abs(logFC) <= log2(fc_soglia) );
    
    mRNAid(ind)=[];
    NmRNA(ind,:)=[];
    CmRNA(ind,:)=[];
    pvalues(ind)=[];
    fdr(ind)=[];
    logFC(ind)=[];
    meanC(ind)=[];
    meanN(ind)=[];

    data=[NmRNA CmRNA];
    pz=[NpzmRNA; CpzmRNA];
    N=length(mRNAid);
    
    filename = strcat(dirnameTxtFilter,'stat-log-dataFiltered.txt');
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', 'gene', 'pvalues', 'fdr', ...
        'mean cancer', 'mean normal', 'log-ratio');
    for i = 1 : N
        fprintf(fid, '%s\t%f\t%f\t%f\t%f\t%f\n', mRNAid{i}, pvalues(i), fdr(i), ...
            meanC(i), meanN(i), logFC(i));
    end
    fclose(fid);
    
    filename = strcat(dirnameTxtFilter,'log-dataFiltered.txt');
    fid = fopen(filename, 'w');
    fprintf(fid, 'gene');
    for c = 1 : numel(pz)
        fprintf(fid, '\t%s', char(pz(c)));
    end
    fprintf(fid, '\n');
    for i = 1 : N
        fprintf(fid, '%s', mRNAid{i});
        for c = 1 : numel(pz)
            fprintf(fid, '\t%f', data(i,c));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
else
    disp('------------------------------------------')
    disp('ATTENZIONE: non hai fatto il log2 dei dati')
    disp('------------------------------------------')
    
    FC = ( meanC + 1) ./ ( meanN + 1);
    
    figure
    [freq,xcenter]=hist(FC,100);
    h=bar(xcenter,freq,'hist');
    
    x1 = abs(xcenter-fc_soglia) < diff(xcenter(1:2))/2; % Find the index of the bin containing the threshold
    x2 = abs(xcenter-(1/fc_soglia)) < diff(xcenter(1:2))/2; % Find the index of the bin containing the threshold
    
    index=(xcenter>=xcenter(x1==1) | xcenter<=xcenter(x2==1));
    colors = [(index(:)+0.5.*(~index(:))) 0.5.*(~index(:)) 0.5.*(~index(:))];
    set(h,'FaceVertexCData',colors);
    xlabel('Ratio',...
        'fontsize',12,'fontname','Times New Roman', 'FontWeight','bold')
    ylabel('Frequency','fontsize',12,'fontname',...
        'Times New Roman', 'FontWeight','bold')
    title('RNA','fontsize',12,'fontname',...
        'Times New Roman', 'FontWeight','bold')
    xlim([min(FC) max(FC)])
    
    ind = find( (FC <= fc_soglia) & (FC >= (1/fc_soglia)) );
    
    mRNAid(ind)=[];
    NmRNA(ind,:)=[];
    CmRNA(ind,:)=[];
    pvalues(ind)=[];
    fdr(ind)=[];
    FC(ind)=[];
    meanC(ind)=[];
    meanN(ind)=[];
    
    data=[NmRNA CmRNA];
    pz=[NpzmRNA; CpzmRNA];
    N=length(mRNAid);
    
    filename = strcat(dirnameTxtFilter,'stat-dataFiltered.txt');
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', 'gene', 'pvalues', 'fdr', ...
        'mean cancer', 'mean normal', 'ratio');
    for i = 1 : N
        fprintf(fid, '%s\t%f\t%f\t%f\t%f\t%f\n', mRNAid{i}, pvalues(i), fdr(i), ...
            meanC(i), meanN(i), FC(i));
    end
    fclose(fid);
    
    filename = strcat(dirnameTxtFilter,'dataFiltered.txt');
    fid = fopen(filename, 'w');
    fprintf(fid, 'gene');
    for c = 1 : numel(pz)
        fprintf(fid, '\t%s', char(pz(c)));
    end
    fprintf(fid, '\n');
    for i = 1 : N
        fprintf(fid, '%s', mRNAid{i});
        for c = 1 : numel(pz)
            fprintf(fid, '\t%f', data(i,c));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end

end