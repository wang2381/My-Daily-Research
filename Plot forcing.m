




%******************************************1000 mean
subplot(3,2,1);
load TS.mat;
clear Ts1000;
load TS_RCP26;
load TS_RCP85;
%%
for i=1:22

Ts1000(i) = mean(Tsmean(i*1000*12-11999:i*1000*12));

end

Ts1000(1,23)=mean(TS_RCP26);
Ts1000_r(1,1:22)=Ts1000(1,1:22);
Ts1000_r(1,23)=mean(TS_RCP85);
Ts1000_c(1,1:23)=Ts1000(1,1:23);
Ts1000_c(1,23)=Ts1000_c(1,23)+1.1;
%Ts1000(23) = mean(Tsmean(23*1000*12-11999:end));

i=11:23;
plot(i,Ts1000(1,i),'-o',i,Ts1000_r(1,i),'-o',i,Ts1000_c(1,i),'-o','linewidth',2,'markersize',10)

x = [11 12 13 14 15 16 17 18 19 20 21 22 23];
x1 = {'12-11','11-10','10-9','9-8','8-7','7-6','6-5','5-4','4-3','3-2','2-1','1-0','21th'};

set(gca,'xtick',x,'fontname','times');
set(gca,'xticklabel',x1,'fontname','times','fontsize',20);
xlabel('Ka','fontsize',20,'fontname','times');
ylabel('Annual Temperature (^{\circ} C)','fontname','times','fontsize',20);
axis([10 24 25 35]);

%*********************************************2000 seasonality
subplot(3,2,2);
clear Ts2000;

for i=1:12

for j=1:11

Ts2000(i,j) = mean(Tsmean([2000*12*(j-1)+i:12:2000*12*j]));

end

end

%for i=1:12

%Ts2000(i,12) = mean(Tsmean([12*2000*11+i:12:end]));

%end




i=1:12

plot(i,Ts2000(:,6),i,Ts2000(:,7),i,Ts2000(:,8),i,Ts2000(:,9),i,Ts2000(:,10),i,Ts2000(:,11),i,TS_RCP26(i),'--',i,TS_RCP26(i)+1.1,'--b', i,TS_RCP85(i),'--r');

legend('12-10k','10-8k','8-6k','6-4k','4-2k','2-0k','RCP 2.6','RCP 4.5','RCP 8.5');
x = [1 2 3 4 5 6 7 8 9 10 11 12];
x1 = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Spt','Oct','Nov','Dec'};

set(gca,'xtick',x,'fontname','times','fontsize',20);
set(gca,'xticklabel',x1,'fontname','times','fontsize',20);
ylabel('Monthly Temperature (^{\circ} C)','fontname','times','fontsize',20);
axis([1 12 22 33]);















%*****************************************************************************Precp (total precipitation)


%******************************************1000 mean


load PP_RCP26;
load PP_RCP85;
load Prec.mat
Precpmean1(120001:120013) = Precpmean(120001:120013);
%Precpmean(1:168000) = Precpmean(1:168000)+10;
%Precpmean(2:12:168000) = Precpmean(2:12:168000)+15;
%Precpmean(6:12:168000) = Precpmean(6:12:168000)-10;
%Precpmean(7:12:168000) = Precpmean(7:12:168000)-10;
%Precpmean(8:12:168000) = Precpmean(8:12:168000)-10;
%Precpmean(9:12:168000) = Precpmean(9:12:168000)-10;
%Precpmean(12:12:168000) = Precpmean(12:12:168000)+15;
%Precpmean(11:12:168000) = Precpmean(11:12:168000)+15;
Precpmean(120001:120013) = Precpmean1(120001:120013);
Precpmean(168007:12:215995) = Precpmean(168007:12:215995)-10;
Precpmean(168008:12:215996) = Precpmean(168008:12:215996)-10;
Precpmean(168009:12:215997) = Precpmean(168009:12:215997)-10;
Precpmean(168010:12:215998) = Precpmean(168010:12:215998)-5;
Precpmean(240007:12:end) = Precpmean(240007:12:end)+5;
Precpmean(240008:12:end) = Precpmean(240008:12:end)+30;
Precpmean(216008:12:239996) = Precpmean(216008:12:239996)+20;
Precpmean(240009:12:end) = Precpmean(240009:12:end)+25;
subplot(3,2,3);
clear Precp1000;

for i=1:22

Precp1000(i) = mean(Precpmean(i*1000*12-11999:i*1000*12))*12;

end
Precp1000(1,23)=sum(PP_RCP85);
Precp1000(2,23)=sum(PP_RCP26);
Precp1000(3,23)=sum(PP_RCP26)+30;
mat=[Precp1000(1,:)',Precp1000(3,:)',Precp1000(2,:)'];
x = [11 12 13 14 15 16 17 18 19 20 21 22 23];
x1 = {'12-11','11-10','10-9','9-8','8-7','7-6','6-5','5-4','4-3','3-2','2-1','1-0','21th'};


bar(mat,'facecolor',[0.7 0.78 1],'edgecolor','k','barwidth',1);



set(gca,'xtick',x,'fontname','times','fontsize',20);
set(gca,'xticklabel',x1,'fontname','times','fontsize',20);
xlabel('Ka','fontname','times','fontsize',20);
ylabel('Annual Precipitation (mm)','fontname','times','fontsize',20);
axis([10 24 2000 3500]);

%*********************************************2000 seasonality
subplot(3,2,4);
clear Precp2000;


for i=1:12

for j=1:11

Precp2000(i,j) = mean(Precpmean([2000*12*(j-1)+i:12:2000*12*j]));

end

end


i=1:12

plot(i,Precp2000(:,6),i,Precp2000(:,7),i,Precp2000(:,8),i,Precp2000(:,9),i,Precp2000(:,10),i,Precp2000(:,11),i,PP_RCP26(i),'--',i,PP_RCP26(i)+6,'--',i,PP_RCP85(i),'--r');

%legend('12-10k','10-8k','8-6k','6-4k','4-2k','2-0k');
x = [1 2 3 4 5 6 7 8 9 10 11 12];
x1 = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Spt','Oct','Nov','Dec'};

set(gca,'xtick',x,'fontname','times');
set(gca,'xticklabel',x1,'fontname','times','fontsize',20);
ylabel('Monthly Precipitation (mm)','fontname','times','fontsize',20);
axis([1 12 50 400])

save Prec1.mat Precpmean






%*********************************************************************************extract par (proportional to TS)
clear solsmean;
load TS_RCP26;
load TS_RCP85;
load TS.mat
solsmean(:,:) = Tsmean(:,:) * 5;

solsmean(1:216000) = solsmean(1:216000)-1;
solsmean(1:168000) = solsmean(1:168000)+2;
solsmean(8:12:168000) = solsmean(8:12:168000)-1;
solsmean(120008) = solsmean(120008)+1;
solsmean(120001:120013) = solsmean(120001:120013)+1;
solsmean(168008:12:215996) = solsmean(168008:12:215996)-1;
solsmean(168009:12:215997) = solsmean(168009:12:215997)-1;
solsmean(168010:12:215998) = solsmean(168010:12:215998)-1;
solsmean(168011:12:215999) = solsmean(168011:12:215999)-1.5;
solsmean(168012:12:216000) = solsmean(168012:12:216000)-1;
solsmean(240008:12:end) = solsmean(240008:12:end)+0.5;
solsmean(216008:12:239996) = solsmean(216008:12:239996)+0.5;
solsmean(240009:12:end) = solsmean(240009:12:end)+1.5;
solsmean(216009:12:239997) = solsmean(216009:12:239997)+1;
solsmean(240010:12:end) = solsmean(240010:12:end)+2.5;
solsmean(240011:12:end) = solsmean(240011:12:end)+1.5;
solsmean(240012:12:end) = solsmean(240012:12:end)+1.5;
solsmean(216010:12:239998) = solsmean(216010:12:239998)+1;

save sols.mat solsmean
%********************************************************************slight revision
solsmean(:,:) = Tsmean(:,:) * 5;


%solsmean(240008:12:end) = solsmean(240008:12:end)+0.5;
%solsmean(216008:12:239996) = solsmean(216008:12:239996)+0.5;
%solsmean(240009:12:end) = solsmean(240009:12:end)+1.5;
%solsmean(216009:12:239997) = solsmean(216009:12:239997)+1;
%solsmean(240010:12:end) = solsmean(240010:12:end)+2.5;
%solsmean(240011:12:end) = solsmean(240011:12:end)+1.5;
%solsmean(240012:12:end) = solsmean(240012:12:end)+1.5;
%solsmean(216010:12:239998) = solsmean(216010:12:239998)+1;





%*****************************************************************************Par (proportional to TS)


%******************************************1000 mean
subplot(3,2,5);

clear sols1000;

for i=1:22

sols1000(i) = mean(solsmean(i*1000*12-11999:i*1000*12));

end

%Ts1000(23) = mean(Tsmean(23*1000*12-11999:end));
sols1000(1,23)=sols1000(1,22)+(mean(TS_RCP85)*5-sols1000(1,22))/2;
sols1000(2,23)=sols1000(1,22)+(mean(TS_RCP26)*5-sols1000(1,22))/2;
sols1000_c(1,:)=sols1000(2,:);
sols1000_c(1,23)=sols1000(2,23)+2;

mat=[sols1000(1,:)',sols1000_c(1,:)',sols1000(2,:)'];

bar(mat,'facecolor',[0.86 0.86 0.86],'edgecolor','k','barwidth',1);

x = [11 12 13 14 15 16 17 18 19 20 21 22 23];
x1 = {'12-11','11-10','10-9','9-8','8-7','7-6','6-5','5-4','4-3','3-2','2-1','1-0','21th'};

set(gca,'xtick',x,'fontname','times');
set(gca,'xticklabel',x1,'fontname','times','fontsize',20);
xlabel('Ka','fontsize',20,'fontname','times');
ylabel('Annual PAR (W m^{-2})','fontname','times','fontsize',20);
axis([10 24 100 150])

%*********************************************2000 seasonality
subplot(3,2,6);
clear sols2000;

for i=1:12

for j=1:11

sols2000(i,j) = mean(solsmean([2000*12*(j-1)+i:12:2000*12*j]));

end

end

%for i=1:12

%Ts2000(i,12) = mean(Tsmean([12*2000*11+i:12:end]));

%end




i=1:12

plot(i,sols2000(:,6),i,sols2000(:,7),i,sols2000(:,8),i,sols2000(:,9),i,sols2000(:,10),i,sols2000(:,11),i,TS_RCP26(i)*5,'--',i,TS_RCP26(i)*5+2,'--',i,TS_RCP85(i)*5-8,'--r');

%legend('12-10k','10-8k','8-6k','6-4k','4-2k','2-0k');
x = [1 2 3 4 5 6 7 8 9 10 11 12];
x1 = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Spt','Oct','Nov','Dec'};

set(gca,'xtick',x,'fontname','times','fontsize',20);
set(gca,'xticklabel',x1,'fontname','times','fontsize',20);
ylabel('Monthly PAR (W m^{-2})','fontname','times','fontsize',20);
axis([1 12 120 160]);