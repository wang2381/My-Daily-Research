%This part does the manipulation of the extracted NC data from the Cluster
%NC files

%input the locations for NC peatlands:

clc
clear
load('C:\research\paper_conference\data_matlab\loc.txt');
clim_ts(1:2356,1:156480)=0;
clim_pp(1:2356,1:156480)=0;
clim_vp(1:2356,1:156480)=0;
clim_par(1:2356,1:156480)=0;

% call functions to obtatin the climate data for each pixel at half by half
% degree
load('C:\research\paper_conference\data_matlab\temp.mat');
load('C:\research\paper_conference\data_matlab\Precipitation.mat');
load('C:\research\paper_conference\data_matlab\vp.mat');
load('C:\research\paper_conference\data_matlab\Par.mat');
%%
for ii=1:2356
    
    [i,j]=findloc(loc(ii,1),loc(ii,2));
    [climate_ts, climate_pp, climate_vp, climate_par]=truncate(i,j,Ts_1,Prec,Vp,Par);
    clm_ts(ii,1:156480)=climate_ts;
    clm_pp(ii,1:156480)=climate_pp;
    clm_vp(ii,1:156480)=climate_vp;
    clm_par(ii,1:156480)=climate_par;
    ts_ij(ii,1:2)=[i,j];
    
end

%save the results due to high computation
save('clm_ts.mat','clm_ts');
save('clm_pp.mat','clm_pp');
save('clm_vp.mat','clm_vp');
save('clm_par.mat','clm_par');
save('ts_ij.mat','ts_ij');
%%
load('clm_ts.mat');
load('clm_pp.mat');
load('clm_vp.mat');
load('clm_par.mat');
load('ts_ij.mat');
%load the RCP (1948-2099) data and intepolate
path='C:\research\paper_conference\data_matlab\na_ts_RCP26_extracted.txt';
[RCP_ij, mjan, mfeb, mmar, mapr, mmay, mjun, mjul, maug, msep, moct, mnov, mdec, mlon, mlat]=RCP(clm_ts,clm_pp,clm_vp,clm_par, path);
%%
%interpolating nc data using RCP data (for ts and pp)

[RCP_diff ts_itpd]=interpolating(RCP_ij,ts_ij, mjan, mfeb, mmar, mapr, mmay, mjun, mjul, maug, msep, moct, mnov, mdec, mlon, mlat, clm_ts);
save('ts_itpd.mat','ts_itpd');



path='C:\research\paper_conference\data_matlab\na_pp_RCP26_extracted.txt';
[RCP_ij, mjan, mfeb, mmar, mapr, mmay, mjun, mjul, maug, msep, moct, mnov, mdec, mlon, mlat]=RCP(clm_ts,clm_pp,clm_vp,clm_par, path);
[RCP_diff pp_itpd]=interpolating(RCP_ij,ts_ij, mjan, mfeb, mmar, mapr, mmay, mjun, mjul, maug, msep, moct, mnov, mdec, mlon, mlat, clm_pp);
save('pp_itpd.mat','pp_itpd');

%%
%interpolating for vp and par

[par_itpd]=par_interpolating(RCP_ij,ts_ij, clm_par);
save('par_itpd.mat','par_itpd');

[vp_itpd]=par_interpolating(RCP_ij,ts_ij, clm_vp);
save('vp_itpd.mat','vp_itpd');

%%
%concatenating the interpolated forcing with the RCP data. 

load('C:\research\paper_conference\data_matlab\ts_itpd.mat');
load('C:\research\paper_conference\data_matlab\pp_itpd.mat');
load('C:\research\paper_conference\data_matlab\par_itpd.mat');
load('C:\research\paper_conference\data_matlab\vp_itpd.mat');


[ts_cont, pp_cont, par_cont, vp_cont]=concatenating(ts_itpd, pp_itpd, par_itpd, vp_itpd);


%%

save('par_cont.mat','par_cont');
save('pp_cont.mat','pp_cont');
save('vp_cont.mat','vp_cont');
save('ts_cont.mat','ts_cont');

