%***********************************************Simulated C Fluxes****************************************
clear aa aaa rate_20 rate_500
FILEID = fopen('C:\research\paper4\PF\NPP_rate.txt');
AA = textscan(FILEID,'%f','delimiter',',');
fclose(FILEID);
FILEID = fopen('C:\research\paper4\PF\RH_rate.txt');
BB = textscan(FILEID,'%f','delimiter',',');
fclose(FILEID);
FILEID = fopen('C:\research\paper4\PF\VSM1_rate.txt');
DD = textscan(FILEID,'%f','delimiter',',');
fclose(FILEID);
FILEID = fopen('C:\research\paper4\PF\wtd.txt');
CC = textscan(FILEID,'%f %f','delimiter',',');
fclose(FILEID);
%%
npp(:,1) = AA{1};
npp1(:,1) = npp(2001:12000,1);

for ii = 1:10000/20

     npprate_20(ii,1) = mean(npp1(((ii-1)*20+1):ii*20,1))*12;
     
end
%rate_20(301:350,1) = rate_20(301:350,1) - 10;

rh(:,1) = BB{1};
rh1(:,1) = rh(2001:12000,1);

for ii = 1:10000/20

     rhrate_20(ii,1) = mean(rh1(((ii-1)*20+1):ii*20,1))*12;
     
end


vsm1(:,1) = DD{1};
vsm11(:,1) = vsm1(2001:12000,1);

for ii = 1:10000/20

     vsm1rate_20(ii,1) = mean(vsm11(((ii-1)*20+1):ii*20,1));
     
end
%%
clear wtd1 wtd2 wtdrate_20
wtd(:,1) = CC{1};
wtd1(:,1) = wtd(2000*12+1:12000*12,1);
for ii = 1:120000/12
wtd2(ii,1) = mean(wtd1((ii-1)+1:ii*12,1));
end

for ii = 1:10000/20

     wtdrate_20(ii,1) = mean(wtd2(((ii-1)*20+1):ii*20,1))-20;
     
end


%%
%**********************Plot*******************************************************
npprate_20_u(:,:)=npprate_20(:,:)+21;
npprate_20_l(:,:)=npprate_20(:,:)-21;
rhrate_20_u(:,:)=rhrate_20(:,:)+21;
rhrate_20_l(:,:)=rhrate_20(:,:)-21;
npprate_20_u(1,:)=nan;
npprate_20_l(1,:)=nan;
rhrate_20_u(1,:)=nan;
rhrate_20_l(2,:)=nan;

subplot(2,2,1);
ii=1:500;
hold on
plot(ii,npprate_20_u(ii,1),'color',[0.95 0.87 0.73],'linewidth',2);
plot(ii,npprate_20_l(ii,1),'color',[0.95 0.87 0.73],'linewidth',2);
plot(ii,npprate_20(ii,1),'color',[0.95 0.87 0.73],'linewidth',2);
x = [1 50 100 150 200 250 300 350 400 450 500];
x1 = {'10','9','8','7','6','5','4','3','2','1','0'};
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',x1,'fontsize',20);
xlabel('Ka','fontsize',20);
ylabel('NPP (g C m^{-2} yr^{-1})','fontsize',20);
axis([1 500 700 850]);
hold off

subplot(2,2,2);
hold on
ii=1:500;
plot(ii,rhrate_20_u(ii,1),'color',[0.73 0.83 0.96],'linewidth',2);
plot(ii,rhrate_20_l(ii,1),'color',[0.73 0.83 0.96],'linewidth',2);
plot(ii,rhrate_20(ii,1),'color',[0.73 0.83 0.96],'linewidth',2);
x = [1 50 100 150 200 250 300 350 400 450 500];
x1 = {'10','9','8','7','6','5','4','3','2','1','0'};
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',x1,'fontsize',20);
xlabel('Ka','fontsize',20);
ylabel('RH (g C m^{-2} yr^{-1})','fontsize',20);
axis([1 500 650 780]);
hold off

subplot(2,2,3);
ii=1:500;
plot(ii,vsm1rate_20(ii,1),'color',[0.23 0.44 0.34],'linewidth',2);
x = [1 50 100 150 200 250 300 350 400 450 500];
x1 = {'10','9','8','7','6','5','4','3','2','1','0'};
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',x1,'fontsize',20);
xlabel('Ka','fontsize',20);
ylabel('VSM (%)','fontsize',20);
axis([1 500 47.5 49]);
hold off

subplot(2,2,4);
ii=1:500;
plot(ii,wtdrate_20(ii,1),'color','k','linewidth',2);
x = [1 50 100 150 200 250 300 350 400 450 500];
x1 = {'10','9','8','7','6','5','4','3','2','1','0'};
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',x1,'fontsize',20);
set(gca,'ydir','reverse');
xlabel('Ka','fontsize',20);
ylabel('Water Table (mm)','fontsize',20);
axis([1 500]);
hold off




