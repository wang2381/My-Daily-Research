% The function extacts the RCP data and interpolates

function [RCP_ij,mjan, mfeb, mmar, mapr, mmay, mjun, mjul, maug, msep, moct, mnov, mdec, mlon, mlat]=RCP(clm_ts,clm_pp,clm_vp,clm_par,path)
clear jan feb mar apr may jun jul aug sep oct nov dec mlon mlat;

FILEID = fopen(path);
A = textscan(FILEID,'%f %f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s','delimiter',',');
fclose(FILEID);


lon(1:355528,1) = A{1};
lat(1:355528,1) = A{2};

jan(1:355528,1) = A{10};
feb(1:355528,1) = A{11};
mar(1:355528,1) = A{12};
apr(1:355528,1) = A{13};
may(1:355528,1) = A{14};
jun(1:355528,1) = A{15};
jul(1:355528,1) = A{16};
aug(1:355528,1) = A{17};
sep(1:355528,1) = A{18};
oct(1:355528,1) = A{19};
nov(1:355528,1) = A{20};
dec(1:355528,1) = A{21};

%do monthly average

for ii=1:2339
    mlon(ii,1) = mean(lon(ii*152-151:ii*152,1));
    mlat(ii,1) = mean(lat(ii*152-151:ii*152,1));
    mjan(ii,1) = mean(jan(ii*152-151:ii*152,1));
    mfeb(ii,1) = mean(feb(ii*152-151:ii*152,1));
    mmar(ii,1) = mean(mar(ii*152-151:ii*152,1));
    mapr(ii,1) = mean(apr(ii*152-151:ii*152,1));
    mmay(ii,1) = mean(may(ii*152-151:ii*152,1));
    mjun(ii,1) = mean(jun(ii*152-151:ii*152,1));
    mjul(ii,1) = mean(jul(ii*152-151:ii*152,1));
    maug(ii,1) = mean(aug(ii*152-151:ii*152,1));
    msep(ii,1) = mean(sep(ii*152-151:ii*152,1));
    moct(ii,1) = mean(oct(ii*152-151:ii*152,1));
    mnov(ii,1) = mean(nov(ii*152-151:ii*152,1));
    mdec(ii,1) = mean(dec(ii*152-151:ii*152,1));
    
    [i,j]=findloc(mlon(ii),mlat(ii));
    RCP_ij(ii,1)=i;
    RCP_ij(ii,2)=j;
end

