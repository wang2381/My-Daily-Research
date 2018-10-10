%temperature********************************************************************************

load('/scratch/conte/w/wang2381/paper4/matlab/TS.mat');
FILEID = fopen('/scratch/conte/w/wang2381/paper4/matlab/pastaza_ts.txt','w');

lon = -75;
lat = -5;
jj=1;
%Tsmean(1,264480-48000:264480) = Tsmean(1,264480-48000:264480) + 2;
for ii=1:12:size(Tsmean,2)-10000*12

C1=zeros(1,12);
C1(1,:)=Tsmean(10000*12+jj*12-11:10000*12+jj*12);
C1_s=sum(C1);
C1_max=max(C1);
C1_m=mean(C1);
C1_min=min(C1);
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f, %s\n',lon,lat, 'TS' ,1439,jj-1,C1_s,C1_max,C1_m,C1_min,C1, 'P');
jj = jj+1;

end
fclose(FILEID);




%Precipitation*********************************************************************************

load('/scratch/conte/w/wang2381/paper4/matlab/Prec1.mat');
FILEID = fopen('/scratch/conte/w/wang2381/paper4/matlab/pastaza_pp.txt','w');
 
%Precpmean(1,264480-48000:264480) = Precpmean(1,264480-48000:264480) + 100;

lon = -75;
lat = -5;
jj=1;

for ii=1:12:size(Precpmean,2)-10000*12

C1=zeros(1,12);
C1(1,:)=Precpmean(10000*12+jj*12-11:10000*12+jj*12);
C1_s=sum(C1);
C1_max=max(C1);
C1_m=mean(C1);
C1_min=min(C1);
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f, %s\n',lon,lat, 'PREC' ,1439,jj-1,C1_s,C1_max,C1_m,C1_min,C1, 'P');
jj = jj+1;

end
fclose(FILEID);




%Vapor pressure in kPa*********************************************************************************

load('/scratch/conte/w/wang2381/paper4/matlab/vp.mat');
FILEID = fopen('/scratch/conte/w/wang2381/paper4/matlab/vpr.txt','w');

lon = -75;
lat = -5;
jj=1;

for ii=1:12:size(vpmean,1)

C1=zeros(1,12);
C1(1,:)=vpmean(jj*12-11:jj*12);
C1_s=sum(C1);
C1_max=max(C1);
C1_m=mean(C1);
C1_min=min(C1);
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f, %s\n',lon,lat, 'VP' ,1439,jj-1,C1_s,C1_max,C1_m,C1_min,C1, 'P');
jj = jj+1;

end
fclose(FILEID);






%incoming shortwave solar radiation (PAR/sols)********************************************************************************

load('/scratch/conte/w/wang2381/paper4/matlab/sols.mat');
FILEID = fopen('/scratch/conte/w/wang2381/paper4/matlab/pastaza_par.txt','w');
%solsmean(1,264480-48000:264480) = solsmean(1,264480-48000:264480) + 2;

lon = -75;
lat = -5;
jj=1;

for ii=1:12:size(solsmean,2)-10000*12

C1=zeros(1,12);
C1(1,:)=solsmean(10000*12+jj*12-11:10000*12+jj*12);
C1_s=sum(C1);
C1_max=max(C1);
C1_m=mean(C1);
C1_min=min(C1);
fprintf(FILEID,'%.1f,%.1f, %s ,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f, %s\n',lon,lat, 'PAR' ,1439,jj-1,C1_s,C1_max,C1_m,C1_min,C1, 'P');
jj = jj+1;

end
fclose(FILEID);

