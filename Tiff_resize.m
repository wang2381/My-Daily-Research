clc
clear
[A,map] = imread('C:\research\paper4\map from draper\Draper_clip_for_elimination','tiff');



B=double(A);
A(A==30)=251;

[C,newmap]=imresize(B,map,0.05,'nearest');
figure;
imshow(A,map);
C(C==7)=0;
%here plot the original clipped map for figure 1 in the paper
%%

for i=1:206
    for j=1:283
        if(B(i,j)==5)
            B(i,j)=0;
        end
        if((j<=280)&&(j>=223)&&(i>=18)&&(i<=55)&&B(i,j)==1)
            B(i,j)=5;
        end
        if((j<=274)&&(j>=240)&&(i>=58)&&(i<=83)&&B(i,j)==1)
            B(i,j)=5;
        end
       if((j<=248)&&(j>=118)&&(i>=25)&&(i<=173)&&B(i,j)==1)
            B(i,j)=5;
        end
    end
end
%%
%revised 5/29 to eliminate the non-grey area in Draper 2014, we output the location outside the PMFB. 
D=C;
D(D(:,:)==7)=0;


%D(D(:,:)==4)=13;
%D(D(:,:)==8)=13;
%D(D(:,:)==15)=13;
FILEID = fopen('C:\research\paper4\map from Draper\loc_seasonal','w');
for i=1:206
    for j=1:283
        if(D(i,j)==0)
        lon = i;
        lat = j;
fprintf(FILEID,'%d,%d\n',lon,lat);
        end
        end
end
fclose(FILEID);



save pastaza_veg.mat D;


%%
figure;
load pastaza_veg.mat;
C = double(B);
C(C(:,:)==0) = 9.46;
C(C(:,:)==3) = 40;
C(C(:,:)==5) = 190;
C(C(:,:)==6) = 110;
C(C(:,:)==4) = nan;


cd C:\research\paper_conference\m_map
intlon=4.5/282;
intlat=3.5/205;
lon=[-77.8:intlon:-73.3];
lat=[-3.45:-intlat:-6.95];
[flon,flat]=meshgrid(lon,lat);
lonmin=-77.8;lonmax=-73.3;
latmin=-6.95;latmax=-3.45;
m_proj('lambert','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_pcolor(flon,flat,C);
m_grid('box','fancy','tickdir','in','fontsize',12);
shading flat;
cd C:\research\paper4\matlab
%*******************************output location for open peatland for
%further use*********************************************************


FILEID = fopen('C:\research\paper4\matlab\loc_OP.txt','w');
kk=1;
for ii=1:206
    for jj=1:283
        if(C(ii,jj)==54)
            loc_OP(kk,1)=ii;
            loc_OP(kk,2)=jj;
            
            fprintf(FILEID,'%.1f,%.1f\n',loc_OP(kk,1),loc_OP(kk,2));
            kk=kk+1;
        end
    end
end
fclose(FILEID);

%%
%output data
%*********************************Palm Swamp******************************
load pastaza_veg.mat;
FILEID = fopen('C:\research\paper4\pastaza_veg_palm.txt','w');
for i=1:206
    for j=1:283
        if(B(i,j)==6||B(i,j)==4)
        lon = i;
        lat = j;
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d, %s , %s\n',lon,lat, 'TEMVEG' ,-111,9, 'Pastaza' , 'Pastaza');
        end
    end
end
fclose(FILEID);

FILEID = fopen('C:\research\paper4\pastaza_elev_palm.txt','w');
for i=1:206
    for j=1:283
        if(B(i,j)==6||B(i,j)==4)
        lon = i;
        lat = j;
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d, %s\n',lon,lat, 'ELEV' ,1439,30, 'Pastaza');
        end
        end
end
fclose(FILEID);


FILEID = fopen('C:\research\paper4\pastaza_soil_palm.txt','w');
for i=1:206
    for j=1:283
        if(B(i,j)==6||B(i,j)==4)
        lon = i;
        lat = j;
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d,%d,%d,%d, %s , %s\n',lon,lat, 'TEXTURE' ,1439,45,40,15,2, 'FAO' , 'Pastaza');
        end
        end
end
fclose(FILEID);




%*********************************Pole Forest******************************


FILEID = fopen('C:\research\paper4\pastaza_veg_forest.txt','w');
for i=1:206
    for j=1:283
        if(B(i,j)==5)
        lon = i;
        lat = j;
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d, %s , %s\n',lon,lat, 'TEMVEG' ,-111,9, 'Pastaza' , 'Pastaza');
        end
    end
end
fclose(FILEID);

FILEID = fopen('C:\research\paper4\pastaza_elev_forest.txt','w');
for i=1:206
    for j=1:283
        if(B(i,j)==5)
        lon = i;
        lat = j;
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d, %s\n',lon,lat, 'ELEV' ,1439,30, 'Pastaza');
        end
        end
end
fclose(FILEID);


FILEID = fopen('C:\research\paper4\pastaza_soil_forest.txt','w');
for i=1:206
    for j=1:283
        if(B(i,j)==5)
        lon = i;
        lat = j;
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d,%d,%d,%d, %s , %s\n',lon,lat, 'TEXTURE' ,1439,45,40,15,2, 'FAO' , 'Pastaza');
        end
        end
end
fclose(FILEID);




%*********************************Terra Firme******************************
load pastaza_veg.mat;
FILEID = fopen('C:\research\paper4\pastaza_veg_up.txt','w');
for i=1:206
    for j=1:283
        if(B(i,j)==3||B(i,j)==0)
        lon = i;
        lat = j;
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d, %s , %s\n',lon,lat, 'TEMVEG' ,-111,9, 'Pastaza' , 'Pastaza');
        end
    end
end
fclose(FILEID);

FILEID = fopen('C:\research\paper4\pastaza_elev_up.txt','w');
for i=1:206
    for j=1:283
        if(B(i,j)==3||B(i,j)==0)
        lon = i;
        lat = j;
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d, %s\n',lon,lat, 'ELEV' ,1439,30, 'Pastaza');
        end
        end
end
fclose(FILEID);


FILEID = fopen('C:\research\paper4\pastaza_soil_up.txt','w');
for i=1:206
    for j=1:283
        if(B(i,j)==3||B(i,j)==0)
        lon = i;
        lat = j;
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d,%d,%d,%d, %s , %s\n',lon,lat, 'TEXTURE' ,1439,45,40,15,2, 'FAO' , 'Pastaza');
        end
        end
end
fclose(FILEID);



%%


