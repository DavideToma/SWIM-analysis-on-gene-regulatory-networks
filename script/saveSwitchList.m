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

function []=saveSwitchList(A,dirnameAllSwitch,switc,mir)

if(mir)
    
    found=strmatch('hsa-',switc);
    switch_mir=switc(found);
    switc(found)=[];
    
    switch_rna={};
    for i=1:numel(switc)
        mat=regexp(switc{i},'\|','split');
        switch_rna{i}=mat{1};
    end
    
    clear switc
    
    switc=[switch_rna'; switch_mir'];
    
    filename=strcat('C:\Users\david\Desktop\SWIM\project\exam_project\all-switch\','lusc','.txt');
    fid=fopen(filename,'w');
    for i=1:numel(switc)
        fprintf(fid,'%s\n',char(switc(i)));
    end
    fclose(fid);
    
else
    
    switch_rna={};
    for i=1:numel(switc)
        mat=regexp(switc{i},'\|','split');
        if(~isempty(mat{1}))
            switch_rna{i}=mat{1};
        else
            switch_rna{i}=switc{i};
        end
    end
    
    switc=switch_rna';
        
    filename=strcat('C:\Users\david\Desktop\SWIM\project\exam_project\all-switch\','lusc','.txt');
    fid=fopen(filename,'w');
    for i=1:numel(switc)
        fprintf(fid,'%s\n',char(switc(i)));
    end
    fclose(fid);
    
end

end

