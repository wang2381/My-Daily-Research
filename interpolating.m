% This function does the interpolation and concatenation

function [RCP_diff clm_data_itpd]=interpolating(RCP_ij, ts_ij, mjan, mfeb, mmar, mapr, mmay, mjun, mjul, maug, msep, moct, mnov, mdec, mlon, mlat, clm_data)
load('C:\research\paper_conference\data_matlab\loc.txt');
load('C:\research\paper_conference\data_matlab\nc_lon.mat');
load('C:\research\paper_conference\data_matlab\nc_lat.mat');
RCP_month(1:2339,1:12)=[mjan, mfeb, mmar, mapr, mmay, mjun, mjul, maug, msep, moct, mnov, mdec];
RCP_ij_copy=RCP_ij;
%obtain the average of the specific region with pixels in it (for RCP)
kk=1;
for ii=1:2339
    if RCP_ij_copy(ii,1)~=-1 && RCP_ij_copy(ii,2)~=-1
        
        RCP_month_tot(kk,1:12)=RCP_month(ii,1:12);
        count=1;
        for jj=ii+1:2339
            if RCP_ij_copy(ii,1)==RCP_ij_copy(jj,1) && RCP_ij_copy(ii,2)==RCP_ij_copy(jj,2)
               RCP_month_tot(kk,1:12)=RCP_month(jj,1:12)+RCP_month_tot(kk,1:12);
               RCP_ij_copy(jj,1:2)=-1;
               count=count+1;
            end
        end
        RCP_month_mean(kk,1:12)=RCP_month_tot(kk,1:12)/count;
        RCP_mean_ij(kk,1)=RCP_ij(ii,1);
        RCP_mean_ij(kk,2)=RCP_ij(ii,2);
        kk=kk+1;
    end
end
%obtain the difference between pixels and the regional average (for RCP)

for ii=1:2339
    for jj=1:size(RCP_mean_ij,1)
        if RCP_ij(ii,1)==RCP_mean_ij(jj,1) && RCP_ij(ii,2)==RCP_mean_ij(jj,2)
            RCP_diff(ii,1:12)=RCP_month(ii,1:12)-RCP_month_mean(jj,1:12);
        end
    end
end

%obtain the difference between pixels and the reginal average (using RCP
%variation to represent the NC paleo data variantion

for ii=1:size(clm_data,1)
    for jj=1:size(RCP_ij,1)
        if ts_ij(ii,1)==RCP_ij(jj,1) && ts_ij(ii,2)==RCP_ij(jj,2) && RCP_ij(jj,1)~=-1
            RCP_ij(jj,1)=-1;
            for kk=1:size(clm_data,2)
                if mod(kk,12)==0
                    clm_data_itpd(jj,kk)=clm_data(ii,kk)+RCP_diff(jj,12);
                else
                    remainder=mod(kk,12);
                    clm_data_itpd(jj,kk)=clm_data(ii,kk)+RCP_diff(jj,remainder);
                end
            end
            break;
        end
            
        
        
    end
end

        

