%*****************************map the soc****************************************


%******************************palm swamp*********************************
clear
figure
subplot(3,2,2);
cd C:\research\paper4\matlab
load C_palm.mat;

C_palm(C_palm(:,:)<=50000)=50000;%60000 for uncertainty calculation
C_palm(C_palm(:,:)>=180000)=180000;

%eliminate the area outside the PMFB
FILENAME = sprintf('C:/research/paper4/map from Draper/loc_seasonal');
FILEID=fopen(FILENAME);
A=textscan(FILEID,'%d%d','delimiter',',');
fclose(FILEID);
loc_lon=A{1};
loc_lat=A{2};
for ii=1:size(loc_lon)
    C_palm(loc_lon(ii),loc_lat(ii))=nan;
end

cd C:\research\Hydrology\Region\m_map
intlon=4.5/282;
intlat=3.5/205;
lon=[-77.6:intlon:-73.1];
lat=[-3.25:-intlat:-6.75];
[flon,flat]=meshgrid(lon,lat);
lonmin=-77.6;lonmax=-73.1;
latmin=-6.75;latmax=-3.25;
m_proj('lambert','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_pcolor(flon,flat,C_palm/1000);
m_grid('box','fancy','tickdir','in','fontsize',12);
shading flat;

set(gca,'clim',[0,250]);
set(gca,'ytick',0:50:250,'yticklabel',{'0','50','100','150','200','250'});
cd C:\research\paper4\matlab
mean_palm=nanmean(nanmean(C_palm(:,:)));
C_palm_tot=8511*1.69*1.69*1000000*mean_palm/1e15;

%***************************open peatland*********************************
subplot(3,2,3);
load C_OP.mat;
C_OP(C_OP(:,:)<=50000*0.49)=50000*0.49;%60000 for uncertainty calculation
C_OP(C_OP(:,:)>=180000*0.49)=180000*0.49;

%eliminate the area outside the PMFB
FILENAME = sprintf('C:/research/paper4/map from Draper/loc_seasonal');
FILEID=fopen(FILENAME);
A=textscan(FILEID,'%d%d','delimiter',',');
fclose(FILEID);
loc_lon=A{1};
loc_lat=A{2};
for ii=1:size(loc_lon)
    C_OP(loc_lon(ii),loc_lat(ii))=nan;
end

cd C:\research\Hydrology\Region\m_map
intlon=4.5/282;
intlat=3.5/205;
lon=[-77.6:intlon:-73.1];
lat=[-3.25:-intlat:-6.75];
[flon,flat]=meshgrid(lon,lat);
lonmin=-77.6;lonmax=-73.1;
latmin=-6.75;latmax=-3.25;
m_proj('lambert','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_pcolor(flon,flat,C_OP/1000);
m_grid('box','fancy','tickdir','in','fontsize',12);
shading flat;

set(gca,'clim',[0,250]);
%set(h,'ytick',0:30:120,'yticklabel',{'0','30','60','90','120'});
cd C:\research\paper4\matlab
mean_open=nanmean(nanmean(C_OP(:,:)));
C_OP_tot=1303*1.69*1.69*1000000*mean_open/1e15;

%******************************pole forest*********************************
clear
subplot(3,2,4);
load C_forest.mat;
C_forest(C_forest(:,:)<=60000)=60000;%100000 for uncertainty calculation
C_forest(C_forest(:,:)>=240000)=240000;

%eliminate the area outside the PMFB
FILENAME = sprintf('C:/research/paper4/map from Draper/loc_seasonal');
FILEID=fopen(FILENAME);
A=textscan(FILEID,'%d%d','delimiter',',');
fclose(FILEID);
loc_lon=A{1};
loc_lat=A{2};
for ii=1:size(loc_lon)
    C_forest(loc_lon(ii),loc_lat(ii))=nan;
end

cd C:\research\Hydrology\Region\m_map
intlon=4.5/282;
intlat=3.5/205;
lon=[-77.6:intlon:-73.1];
lat=[-3.25:-intlat:-6.75];
[flon,flat]=meshgrid(lon,lat);
lonmin=-77.6;lonmax=-73.1;
latmin=-6.75;latmax=-3.25;
m_proj('lambert','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_pcolor(flon,flat,C_forest/1000);
m_grid('box','fancy','tickdir','in','fontsize',12);
shading flat;

set(gca,'clim',[0,250]);
%set(h,'ytick',0:40:280,'yticklabel',{'0','40','80','120','160','200','240','280'});
cd C:\research\paper4\matlab
mean_forest=nanmean(nanmean(C_forest(:,:)));
C_PF_tot=1022*1.69*1.69*1000000*mean_forest/1e15;


%******************************upland*********************************
clear
subplot(3,2,1);
cd C:\research\paper4\matlab
load C_up.mat;
kk=0;
for ii=238:283
    for jj=1:206
        if(isnan(C_up(jj,ii))==0)
            while(isnan(C_up(jj,238+237-ii-kk))==1)
                kk=kk+1;
            end
            C_up(jj,ii)=C_up(jj,238+237-ii-kk);
            if(C_up(jj,ii)<=8930)
            C_up(jj,ii)=8940;
            end
            kk=0;
        end
                
     end
 end
kk=0;
for ii=10:-1:1
    for jj=1:283
        if(isnan(C_up(ii,jj))==0)
            while(isnan(C_up(21-ii+kk,jj))==1)
                kk=kk+1;
            end
            C_up(ii,jj)=C_up(21-ii+kk,jj);
            kk=0;
        end
                
     end
end
            


%eliminate the area outside the PMFB
FILENAME = sprintf('C:/research/paper4/map from Draper/loc_seasonal');
FILEID=fopen(FILENAME);
A=textscan(FILEID,'%d%d','delimiter',',');
fclose(FILEID);
loc_lon=A{1};
loc_lat=A{2};
for ii=1:size(loc_lon)
    C_up(loc_lon(ii),loc_lat(ii))=nan;
end
C_up(1,1)=1000*250;
cd C:\research\Hydrology\Region\m_map
intlon=4.5/282;
intlat=3.5/205;
lon=[-77.6:intlon:-73.1];
lat=[-3.25:-intlat:-6.75];
[flon,flat]=meshgrid(lon,lat);
lonmin=-77.6;lonmax=-73.1;
latmin=-6.75;latmax=-3.25;
m_proj('lambert','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_pcolor(flon,flat,C_up/1000);
m_grid('box','fancy','tickdir','in','fontsize',12);
shading flat;
T = [
    191, 0, 255
    128, 0, 255
    64, 0, 255
    0, 128, 255
    0, 191, 255
    0, 255, 255
    
    0, 255, 128
    0, 255, 0
    191, 255, 0
    255, 255, 0
    255, 191, 0
    255, 128, 0
    255, 0, 0
    ];
 x = [
     8.8
     8.9
     9
     9.1
     9.2
     9.3
     
     20
     50
     70
     100
     150
     200
     250
];
 x2=8.8:0.1:250;
 map = interp1(x,T,x2);
% map(6:2500,1:3)=nan;
 map1=map/255;
colormap(map1);
h=colorbar;
%set(h,'clim',[8.8,9.3]);
y=[8.8,8.9,9.0,9.1,9.2,9.3 ];
set(h,'ytick',y);

%set(map,'ytick',0:50:250);
%set(h,'ytick',8.8:0.2:250,'yticklabel',{'8.8','9.0','9.2','9.4','9.6'});
cd C:\research\paper4\matlab
mean_up=nanmean(nanmean(C_up(:,:)));
C_up_tot=26937*1.69*1.69*1000000*mean_up/1e15;



%**********************************four figures combination***************
subplot(3,2,5.5)
clear
load C_palm.mat;
load C_OP.mat
load C_forest.mat
load C_up.mat

%eliminate the area outside the PMFB
FILENAME = sprintf('C:/research/paper4/map from Draper/loc_seasonal');
FILEID=fopen(FILENAME);
A=textscan(FILEID,'%d%d','delimiter',',');
fclose(FILEID);
loc_lon=A{1};
loc_lat=A{2};
for ii=1:size(loc_lon)
    C_up(loc_lon(ii),loc_lat(ii))=nan;
end

%eliminate the area outside the PMFB
FILENAME = sprintf('C:/research/paper4/map from Draper/loc_seasonal');
FILEID=fopen(FILENAME);
A=textscan(FILEID,'%d%d','delimiter',',');
fclose(FILEID);
loc_lon=A{1};
loc_lat=A{2};
for ii=1:size(loc_lon)
    C_palm(loc_lon(ii),loc_lat(ii))=nan;
end

%eliminate the area outside the PMFB
FILENAME = sprintf('C:/research/paper4/map from Draper/loc_seasonal');
FILEID=fopen(FILENAME);
A=textscan(FILEID,'%d%d','delimiter',',');
fclose(FILEID);
loc_lon=A{1};
loc_lat=A{2};
for ii=1:size(loc_lon)
    C_forest(loc_lon(ii),loc_lat(ii))=nan;
end

%eliminate the area outside the PMFB
FILENAME = sprintf('C:/research/paper4/map from Draper/loc_seasonal');
FILEID=fopen(FILENAME);
A=textscan(FILEID,'%d%d','delimiter',',');
fclose(FILEID);
loc_lon=A{1};
loc_lat=A{2};
for ii=1:size(loc_lon)
    C_OP(loc_lon(ii),loc_lat(ii))=nan;
end
C_palm(C_palm(:,:)<=50000)=50000;
C_palm(C_palm(:,:)>=180000)=180000;
C_OP(C_OP(:,:)<=50000*0.49)=50000*0.49;
C_OP(C_OP(:,:)>=180000*0.49)=180000*0.49;
C_forest(C_forest(:,:)<=30000)=30000;
C_forest(C_forest(:,:)>=240000)=240000;

C_palm(isnan(C_palm(:,:))==1)=0;
C_OP(isnan(C_OP(:,:))==1)=0;
C_forest(isnan(C_forest(:,:))==1)=0;
C_up(isnan(C_up(:,:))==1)=0;

C_tot(1:206,1:283)=0;
C_tot(:,:)=C_palm(:,:)+C_OP(:,:)+C_forest(:,:)+C_up(:,:);
C_tot(C_tot(:,:)==0)=nan;

cd C:\research\Hydrology\Region\m_map
intlon=4.5/282;
intlat=3.5/205;
lon=[-77.6:intlon:-73.1];
lat=[-3.25:-intlat:-6.75];
[flon,flat]=meshgrid(lon,lat);
lonmin=-77.6;lonmax=-73.1;
latmin=-6.75;latmax=-3.25;
m_proj('lambert','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_pcolor(flon,flat,C_tot/1000);
m_grid('box','fancy','tickdir','in','fontsize',12);
shading flat;
T = [
    191, 0, 255
    128, 0, 255
    64, 0, 255
    0, 128, 255
    0, 191, 255
    0, 255, 255
    
    0, 255, 128
    0, 255, 0
    191, 255, 0
    255, 255, 0
    255, 191, 0
    255, 128, 0
    255, 0, 0
    ];
 x = [
     8.8
     8.9
     9
     9.1
     9.2
     9.3
     
     20
     50
     70
     100
     150
     200
     250
];
 x2=8.8:0.1:250;
 map = interp1(x,T,x2);
% map(6:2500,1:3)=nan;
 map1=map/255;
colormap(map1);
h=colorbar;
%set(h,'clim',[8.8,9.3]);
y=[9.0, 15, 20, 50, 100, 150, 200, 250];
set(h,'ytick',y);

cd C:\research\paper4\matlab

%%
C_palm(isnan(C_palm(:,:))==1)=0;
C_OP(isnan(C_OP(:,:))==1)=0;
C(:,:)=C_OP+C_palm;
sum(sum(C==0))
