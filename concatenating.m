function [ts_cont, pp_cont, par_cont, vp_cont]=concatenating(ts_itpd, pp_itpd, par_itpd, vp_itpd)
%for ts
path='C:\research\paper_conference\data_matlab\na_ts_RCP26_extracted.txt';
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
RCP_month(1:355528,1:12)=[jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec];

for ii=1:2339
    ts_cont(ii,:)=[ts_itpd(ii,1:156480), reshape(RCP_month((ii-1)*152+44:ii*152, :)',1,[])];
end


%for pp
path='C:\research\paper_conference\data_matlab\na_pp_RCP26_extracted.txt';
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
RCP_month(1:355528,1:12)=[jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec];

for ii=1:2339
    pp_cont(ii,:)=[pp_itpd(ii,1:156480), reshape(RCP_month((ii-1)*152+44:ii*152, :)',1,[])];
end 


%for par
for ii=1:2339
    par_cont(ii,:)=[par_itpd(ii,1:156480), par_itpd(ii,156480-12*109+1:end)];
end 

%for vp
for ii=1:2339
    vp_cont(ii,:)=[vp_itpd(ii,1:156480), vp_itpd(ii,156480-12*109+1:end)];
end