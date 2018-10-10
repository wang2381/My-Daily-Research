%create locations of North America
clc
clear
cd C:\research\paper_conference\data_matlab

FILENAME=sprintf('loc.txt');
FILEID=fopen(FILENAME,'w');
load C:\research\paper_conference\data_matlab\peatmap.mat;

[a b]=size(peatmap);

row=1;
for i=1:b
    for j=a:-1:1
        if isnan(peatmap(j,i))~=1
            loc(row,1)=-180.0+0.5*(i-1);
            loc(row,2)=10.0+0.5*(a-j);
            row=row+1;
        end
    end
end

for ii=1:2356
fprintf(FILEID,'%.1f,%.1f\n',loc(ii,1),loc(ii,2));
end