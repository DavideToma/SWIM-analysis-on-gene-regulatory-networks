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

function [switc,nodeNameG,APCC,APCCswitc]=GA_DateParty(A,rho,idx,adjMatrix,...
    nodeName,RNAname,nodi,k,dirnameMatSwitch,dirnameTxtSwitch,...
    dirnameFigureSwitch)

M=numel(nodeName);
C=zeros(M,M);
for i=1:M
    for j=1:M
        if (idx(i)==idx(j))
            C(i,j)=1;
        end
    end
end
C(logical(eye(M))) = 0;
zc=C.*adjMatrix;

deg=zeros(M,1);
ideg=zeros(M,1);
for i=1:M
    deg(i)=sum(adjMatrix(i,:));
    ideg(i)=sum(zc(i,:));
end

adegC=zeros(k,1);
stdC=zeros(k,1);
for i=1:k
    tmp=find(idx==i);
    adegC(i)=mean(deg(tmp));
    stdC(i)=std(deg(tmp));
end

P=[];
z=[];
idegG=[];
degG=[];
nodiG=[];
nodeNameG={};
count=0;
for i=1:M
    if (deg(i)>0 && ideg(i)>0)
        count=count+1;
        z(count)=(ideg(i)-adegC(idx(i)))/stdC(idx(i));
        P(count)=1-(ideg(i)/deg(i))^2;
        nodiG(count)=nodi(i);
        nodeNameG(count)=nodeName(i);
        idegG(count)=ideg(i);
        degG(count)=deg(i);
    end
end

attribute1={};
attribute2={};
attribute3={};
for i=1:numel(nodiG)
    if (z(i)<2.5)
        attribute1{i}='non local hub';
        if (P(i)<=0.04)
            attribute2{i}='R1';
            attribute3{i}='Ultra-peripheral nodes';
        elseif (P(i)>0.04  && P(i)<=0.625)
            attribute2{i}='R2';
            attribute3{i}='Peripheral nodes';
        elseif (P(i)>0.625 && P(i)<=0.8)
            attribute2{i}='R3';
            attribute3{i}='Non-hub connectors';
        elseif (P(i)>0.8)
            attribute2{i}='R4';
            attribute3{i}='Non-hub kinless nodes';
        end
    elseif (z(i)>=2.5)
        attribute1{i}='local hub';
        if (P(i)<=0.3)
            attribute2{i}='R5';
            attribute3{i}='Provincial hubs';
        elseif (P(i)>0.3  && P(i)<=0.75)
            attribute2{i}='R6';
            attribute3{i}='Connector hubs';
        elseif (P(i)>0.75)
            attribute2{i}='R7';
            attribute3{i}='Kinless hubs';
        end
    else
        attribute1{i}='z-NAN';
        attribute2{i}='z-NAN';
        attribute3{i}='z-NAN';
    end
end

APCC_all=cell(M,3);
for i=1:M
    
    nearest_rete=nodeName(find(adjMatrix(i,:)>0));
    num_nearest=numel(nearest_rete);
    ind_source=strmatch(nodeName(i),RNAname,'exact');
    
    tmp=[];
    ind_target=zeros(num_nearest,1);
    for j=1:num_nearest
        ind_target(j)=strmatch(nearest_rete(j),RNAname,'exact');
        tmp(j)=rho(ind_source,ind_target(j));
    end
    
    APCC_all{i,1}=mean(tmp);
    APCC_all{i,2}=num_nearest;
    
    if ( (mean(tmp)>=0) && (mean(tmp)<0.5) && deg(i)>=5 )
        APCC_all{i,3}='DATE';
    elseif ( (mean(tmp)>=0.5) && deg(i)>=5 )
        APCC_all{i,3}='PARTY';
    elseif ( (mean(tmp)<0) && deg(i)>=5 )
        APCC_all{i,3}='FIGHT CLUB';
    else 
        APCC_all{i,3}='no hub';
    end
    clear tmp
end

APCC=cell(numel(nodiG),3);
meanAPCC=zeros(numel(nodiG),1);
count=0;
for i=1:numel(nodiG)
    found=strmatch(nodeNameG{i},nodeName,'exact');
    if(~isempty(found))
        count=count+1;
        APCC(count,:)=APCC_all(found,:);
        meanAPCC(count)=APCC_all{found,1};
    end
end

save(strcat(dirnameMatSwitch,'attribute'),'nodeNameG','nodeName',...
    'attribute1','attribute2','attribute3','z','P',...
    'degG','deg','APCC_all','APCC','meanAPCC') 

fid = fopen(strcat(dirnameTxtSwitch,'attribute.txt'),'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
    'nodeName','Hub','Region','type','Degree',...
    'APCC','Date-Party','P','z');
for i=1:numel(nodiG)
    fprintf(fid,'%s\t%s\t%s\t%s\t%d\t%f\t%s\t%f\t%f\n',...
      char(nodeNameG(i)),char(attribute1(i)),...
      char(attribute2(i)),char(attribute3(i)),APCC{i,2},...
      APCC{i,1},APCC{i,3},P(i),z(i));
end
fclose(fid);

switc={};
APCCswitc={};
attributeswitc={};
Pswitc=[];
zswitc=[];
count=0;
for i=1:numel(nodiG)
    if ( strcmp(attribute2(i),'R4') && strcmp(APCC{i,3},'FIGHT CLUB') )
        count=count+1;
        switc(count)=nodeNameG(i);
        attributeswitc(count,:) =[attribute1(i) attribute2(i) attribute3(i)];
        APCCswitc(count,:)=APCC(i,:);
        Pswitc(count)=P(i);
        zswitc(count)=z(i);
    end
end

fid=fopen(strcat(dirnameTxtSwitch,'switch.txt'),'w');
for i=1:numel(switc)
    fprintf(fid,'%s\n',char(switc(i)));
end
fclose(fid);

save(strcat(dirnameMatSwitch,'attribute-switch'),'switc',...
    'attributeswitc','zswitc','Pswitc','APCCswitc')

fid = fopen(strcat(dirnameTxtSwitch,'attribute-switch.txt'),'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
    'nodeName','Hub','Region','type','Degree',...
    'APCC','Date-Party','P','z');
for i=1:numel(switc)
    fprintf(fid,'%s\t%s\t%s\t%s\t%d\t%f\t%s\t%f\t%f\n',...
        char(switc(i)),char(attributeswitc(i,1)),...
        char(attributeswitc(i,2)),char(attributeswitc(i,3)),APCCswitc{i,2},...
        APCCswitc{i,1},APCCswitc{i,3},Pswitc(i),zswitc(i));
end
fclose(fid);

figure('Color', 'w')
scatter(P',z',30,meanAPCC,'filled')
colorbar; caxis([-1 1]);
set(gca,'XTick',0:0.2:1,'fontsize',20)
set(gca,'YTick',-8:2:5,'fontsize',20)
title(upper(A),'fontsize',20,'fontweight','b');
xlabel('Clusterphobic coefficient, K_{\pi}','fontsize',20,'fontweight','b')
ylabel('Whithin-module degree, z_g','fontsize',20,'fontweight','b')
hold on
plot([0,1],[2.5,2.5],'k-','LineWidth',2)
plot([0.04,0.04],[-5,2.5],'k-','LineWidth',2)
plot([0.625,0.625],[-5,2.5],'k-','LineWidth',2)
plot([0.8,0.8],[-5,2.5],'k-','LineWidth',2)
plot([0.3,0.3],[2.5,5],'k-','LineWidth',2)
plot([0.75,0.75],[2.5,5],'k-','LineWidth',2)
text(0.02,-4,'\leftarrow R1',...
    'HorizontalAlignment','left','fontsize',20,'fontweight','b')
text(0.3,-4,'R2',...
    'HorizontalAlignment','left','fontsize',20,'fontweight','b')
text(0.7,-4,'R3',...
    'HorizontalAlignment','left','fontsize',20,'fontweight','b')
text(0.9,-4,'R4',...
    'HorizontalAlignment','left','fontsize',20,'fontweight','b')
text(0.06,4.5,'R5',...
    'HorizontalAlignment','left','fontsize',20,'fontweight','b')
text(0.5,4.5,'R6',...
    'HorizontalAlignment','left','fontsize',20,'fontweight','b')
text(0.9,4.5,'R7',...
    'HorizontalAlignment','left','fontsize',20,'fontweight','b')

saveas(gcf,strcat(dirnameFigureSwitch,'GA.fig'))
saveas(gcf,strcat(dirnameFigureSwitch,'GA.png'))


end