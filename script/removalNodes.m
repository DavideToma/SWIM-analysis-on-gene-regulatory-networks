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

function [] = removalNodes(N,APCC,APCCswitc,switc,adjMatrix,nodeName,...
    nodeNameG,dirnameMatSwitch,dirnameFigureSwitch)

[size0,adj0,nodeName0,degree0,lambda0]=connectivity(adjMatrix, nodeName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
fprintf('Evaluating random attack\n')

ind=randperm(size0);

randomName=nodeName0(ind);

[randomTOremove_list] = removalNodes_list(N,nodeName0,randomName,[],[]);

[fR,sizeR,lambdaR] = connectivityChange(N,size0,randomTOremove_list,...
    adj0,nodeName0);

save(strcat(dirnameMatSwitch,'removalRandomNodes'),'fR','sizeR','lambdaR')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Evaluating targeted attack\n')

[val,ind]=sort(degree0,'descend');

hubTOremove_list=ind(1:N);

[fH,sizeH,lambdaH] = connectivityChange(N,size0,hubTOremove_list,...
    adj0,nodeName0);

save(strcat(dirnameMatSwitch,'removalHubs'),'fH','sizeH','lambdaH')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Evaluating non-switch removal \n')

hub_list=nodeName0(ind);
hub_degree=val;

[name,index]=setdiff(hub_list,switc);

[noswitchhubTOremove_list] = removalNodes_list(N,nodeName0,name,...
    [],hub_degree(index));

[fNS,sizeNS,lambdaNS] = connectivityChange(N,size0,...
    noswitchhubTOremove_list,adj0,nodeName0);

save(strcat(dirnameMatSwitch,'removalNonSwitch'),'fNS','sizeNS','lambdaNS')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Evaluating switch removal \n')

degree_switch=cell2mat(APCCswitc(:,2));

[switchTOremove_list] = removalNodes_list(N,nodeName0,switc,[],degree_switch);

[fS,sizeS,lambdaS] = connectivityChange(N,size0,switchTOremove_list,...
    adj0,nodeName0);

save(strcat(dirnameMatSwitch,'removalSwitch'),'fS','sizeS','lambdaS')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Evaluating fight-club removal \n')

ind_fightclub=strmatch('FIGHT CLUB',APCC(:,3),'exact');

degree_fightclub=cell2mat(APCC(ind_fightclub,2));

[fightclubTOremove_list] = removalNodes_list(N,nodeName0,nodeNameG,...
    ind_fightclub,degree_fightclub);

[fFC,sizeFC,lambdaFC] = connectivityChange(N,size0,fightclubTOremove_list,...
    adj0,nodeName0);

save(strcat(dirnameMatSwitch,'removalFightClub'),'fFC','sizeFC','lambdaFC')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Evaluating date removal \n')

ind_date=strmatch('DATE',APCC(:,3),'exact');

if(~isempty(ind_date))
    degree_date=cell2mat(APCC(ind_date,2));
    
    [dateTOremove_list] = removalNodes_list(N,nodeName0,nodeNameG,ind_date,...
        degree_date);
    
    [fD,sizeD,lambdaD] = connectivityChange(N,size0,dateTOremove_list,...
        adj0,nodeName0);
    
    save(strcat(dirnameMatSwitch,'removalDateHubs'),'ind_date',...
        'fD','sizeD','lambdaD')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Evaluating party removal \n')

ind_party=strmatch('PARTY',APCC(:,3),'exact');

if(~isempty(ind_party))
    degree_party=cell2mat(APCC(ind_party,2));
    
    [partyTOremove_list] = removalNodes_list(N,nodeName0,nodeNameG,ind_party,...
        degree_party);
    
    [fP,sizeP,lambdaP] = connectivityChange(N,size0,partyTOremove_list,...
        adj0,nodeName0);    
    
    save(strcat(dirnameMatSwitch,'removalPartyHubs'),'ind_party',...
        'fP','sizeP','lambdaP')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
removalNodes_plot(ind_date,ind_party,fH,lambdaH,...
    fFC,lambdaFC,fD,lambdaD,fP,lambdaP,fR,lambdaR,fNS,lambdaNS,...
    fS,lambdaS,dirnameFigureSwitch)

end









