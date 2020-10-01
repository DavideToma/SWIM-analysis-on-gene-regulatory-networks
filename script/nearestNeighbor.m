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

function []=nearestNeighbor(data,tissue,logFC,CpzmRNA,NpzmRNA,RNAname,...
    nodeName,switc,nodeNameG,APCC,idx,adjMatrix,dirnameMatSwitch,...
    dirnameTxtSwitch,dirnameFigureSwitch)

N=numel(switc);
M=size(data,2);

APCC_switc=cell(N,3);
for i=1:N
    found=strmatch(switc(i),nodeNameG,'exact');
    APCC_switc(i,1)=APCC(found,1);
    APCC_switc(i,2)=APCC(found,2);
    APCC_switc(i,3)=APCC(found,3);
end

corr_st_neg=[];
nn_name_neg_switch={};
cluster_target_neg=[];
corr_st_pos=[];
nn_name_pos_switch={};
cluster_target_pos=[];
cluster_source=zeros(N,1);
media_corr_pos=zeros(N,1);
media_corr_neg=zeros(N,1);
data_switch=zeros(N,M);
delta=zeros(N,2);
for i=1:N
    found=strmatch(switc(i),nodeName,'exact');
    found1=strmatch(switc(i),RNAname,'exact');
    
    nn_pos=find(adjMatrix(found,:)>0);
    
    k=numel(nn_pos);
    nn_name=nodeName(nn_pos);
    
    nn_data=zeros(k,M);
    nn_fc=zeros(k,1);
    count=0;
    cont=0;
    tmp_neg=[];
    tmp_pos=[];
    for j=1:k
        
        found2=strmatch(nn_name(j),RNAname,'exact');
        
        nn_data(j,:)=data(found2,:);
        nn_fc(j)=logFC(found2);
        
        tmp=corr(data(found2,:)',data(found1,:)','rows','pairwise');
        if( tmp < 0 )
            count=count+1;
            nn_name_neg_switch(i,count)=nn_name(j);
            corr_st_neg(i,count)=tmp;
            tmp_neg(count)=tmp;
            cluster_target_neg(i,count)=idx(nn_pos(j));
        elseif(tmp > 0)
            cont=cont+1;
            nn_name_pos_switch(i,cont)=nn_name(j);
            corr_st_pos(i,cont)=tmp;
            tmp_pos(cont)=tmp;
            cluster_target_pos(i,cont)=idx(nn_pos(j));
        end
    end
    
    cluster_source(i)=idx(found);
    data_switch(i,:)=data(found1,:);
    delta(i,1)=logFC(found1);
    delta(i,2)=mean(nn_fc);
    
    media_corr_neg(i)=mean(tmp_neg);
    media_corr_pos(i)=mean(tmp_pos);
    
end

makeClustergram(data_switch,switc,tissue,CpzmRNA,NpzmRNA,...
    dirnameFigureSwitch,'heatmap_switch','switch')

save(strcat(dirnameMatSwitch,'switc'),'switc', 'data_switch',...
    'tissue','CpzmRNA','NpzmRNA')

a={};
c_tmp=[];
sum_nn_neg_switch=zeros(N,1);
s=0;
for i=1:N
    somma_tmp=0;
    for j=1:size(nn_name_neg_switch,2)
        if(~isempty(nn_name_neg_switch{i,j}))
            s=s+1;
            somma_tmp=somma_tmp+1;
            a{s}=nn_name_neg_switch{i,j};
            c_tmp(s)=cluster_target_neg(i,j);
        end
    end
    sum_nn_neg_switch(i)=somma_tmp;
end

unique_list_nn_neg_switch=unique(a);
unique_list_cluster_nn_neg_switch=unique(c_tmp);

a={};
c_tmp=[];
sum_nn_pos_switch=zeros(N,1);
s=0;
for i=1:N
    somma_tmp=0;
    for j=1:size(nn_name_pos_switch,2)
        if(~isempty(nn_name_pos_switch{i,j}))
            s=s+1;
            somma_tmp=somma_tmp+1;
            a{s}=nn_name_pos_switch{i,j};
            c_tmp(s)=cluster_target_pos(i,j);
        end
    end
    sum_nn_pos_switch(i)=somma_tmp;
end

unique_list_nn_pos_switch=unique(a); 
unique_list_cluster_nn_pos_switch=unique(c_tmp);

unique_list_cluster_switch=unique(cluster_source);

save(strcat(dirnameMatSwitch,'nearestNeighbor'),'sum_nn_neg_switch',...
    'unique_list_nn_neg_switch','unique_list_cluster_nn_neg_switch',...
    'unique_list_cluster_switch','sum_nn_pos_switch',...
    'unique_list_nn_pos_switch','unique_list_cluster_nn_pos_switch',...
    'nn_name_neg_switch','nn_name_pos_switch','switc')

filename = strcat(dirnameTxtSwitch,'nn_name_neg_switch.txt');
fid = fopen(filename, 'w');
fprintf(fid, 'switch');
for c = 1 : size(nn_name_neg_switch,2)
    fprintf(fid, '\t%d', c);
end
fprintf(fid, '\n');
for i = 1 : N
    fprintf(fid, '%s', char(switc{i}));
    for c = 1 : size(nn_name_neg_switch,2)
        fprintf(fid, '\t%s', char(nn_name_neg_switch{i,c}));
    end
    fprintf(fid, '\n');
end
fclose(fid);

filename = strcat(dirnameTxtSwitch,'corr_nn_neg_switch.txt');
fid = fopen(filename, 'w');
fprintf(fid, '%s\t%s', 'switch', 'media corr neg');
for c = 1 : size(nn_name_neg_switch,2)
    fprintf(fid, '\t%d', c);
end
fprintf(fid, '\n');
for i = 1 : N
    fprintf(fid, '%s\t%f', char(switc{i}), media_corr_neg(i));
    for c = 1 : size(nn_name_neg_switch,2)
        fprintf(fid, '\t%f', corr_st_neg(i,c));
    end
    fprintf(fid, '\n');
end
fclose(fid);

filename = strcat(dirnameTxtSwitch,'cluster_nn_neg_switch.txt');
fid = fopen(filename, 'w');
fprintf(fid, '%s', 'cluster switch');
for c = 1 : size(nn_name_neg_switch,2)
    fprintf(fid, '\t%d', c);
end
fprintf(fid, '\n');
for i = 1 : N
    fprintf(fid, '%d', cluster_source(i));
    for c = 1 : size(nn_name_neg_switch,2)
        fprintf(fid, '\t%d', cluster_target_neg(i,c));
    end
    fprintf(fid, '\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%

filename = strcat(dirnameTxtSwitch,'fc_switch.txt');
fid = fopen(filename, 'w');
fprintf(fid, '%s\t%s\n', 'switch', 'log-ratio');
for i = 1 : N
    fprintf(fid, '%s\t%f\n', char(switc{i}), delta(i,1));
end
fclose(fid);

filename = strcat(dirnameTxtSwitch,'nn_name_pos_switch.txt');
fid = fopen(filename, 'w');
fprintf(fid, 'switch');
for c = 1 : size(nn_name_pos_switch,2)
    fprintf(fid, '\t%d', c);
end
fprintf(fid, '\n');
for i = 1 : N
    fprintf(fid, '%s', char(switc{i}));
    for c = 1 : size(nn_name_pos_switch,2)
        fprintf(fid, '\t%s', char(nn_name_pos_switch{i,c}));
    end
    fprintf(fid, '\n');
end
fclose(fid);

filename = strcat(dirnameTxtSwitch,'corr_nn_pos_switch.txt');
fid = fopen(filename, 'w');
fprintf(fid, '%s\t%s', 'switch', 'media corr pos');
for c = 1 : size(nn_name_pos_switch,2)
    fprintf(fid, '\t%d', c);
end
fprintf(fid, '\n');
for i = 1 : N
    fprintf(fid, '%s\t%f', char(switc{i}), media_corr_pos(i));
    for c = 1 : size(nn_name_pos_switch,2)
        fprintf(fid, '\t%f', corr_st_pos(i,c));
    end
    fprintf(fid, '\n');
end
fclose(fid);

filename = strcat(dirnameTxtSwitch,'cluster_nn_pos_switch.txt');
fid = fopen(filename, 'w');
fprintf(fid, '%s', 'cluster switch');
for c = 1 : size(nn_name_pos_switch,2)
    fprintf(fid, '\t%d', c);
end
fprintf(fid, '\n');
for i = 1 : N
    fprintf(fid, '%d', cluster_source(i));
    for c = 1 : size(nn_name_pos_switch,2)
        fprintf(fid, '\t%d', cluster_target_pos(i,c));
    end
    fprintf(fid, '\n');
end
fclose(fid);

end




