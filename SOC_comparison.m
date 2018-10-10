

clear
%***************************aaa observation Peru************************************************************
clear Q Q_10000 Q_20 Q_500 R R_10000 R_20 R_500 S S_10000 S_20 S_500 B B_10000 B_20 B_500 C C_10000 C_20 C_500 A A_10000 A_20 A_500;

%Quistococha
Q(1,1:165) = nan; 
Q(1,165+(1:230)) = 45;
Q(1,165+(231:230+258)) = 125;
Q(1,165+(231+258:231+258+857)) = 52;
Q(1,165+(232+258+857:2197)) = 74;
Q(1,165+(2197:2335)) = nan;
Q_10000(1,1:7500) = nan;
Q_10000(1,7501:10000) = Q(1,:); 
for ii = 1:10000/20
Q_20(1,ii) = nanmean(Q_10000(1,(ii-1)*20+1:ii*20));
end
for ii = 1:10000/500
Q_500(1,ii) = nanmean(Q_10000(1,(ii-1)*500+1:ii*500));
Q_500std(1,ii) = nanstd(Q_10000(1,(ii-1)*500+1:ii*500));
end


%Rinon
R(1,1:385) = nan;
R(1,385+(1:393)) = 126;
R(1,385+(394:849)) = 44;
R(1,385+(850:1392)) = 26;
R(1,385+(1393:1615)) = nan;
R_10000(1,1:8000) = nan;

R_10000(1,8001:10000) = R(1,:);
for ii = 1:10000/20
R_20(1,ii) = nanmean(R_10000(1,(ii-1)*20+1:ii*20));
end
for ii = 1:10000/500
R_500(1,ii) = nanmean(R_10000(1,(ii-1)*500+1:ii*500));
R_500std(1,ii) = nanstd(R_10000(1,(ii-1)*500+1:ii*500));
end

%San Jorge
S(1,1:55) = nan;
S(1,55+(1:205)) = 195;
S(1,55+(206:609)) = 128;
S(1,55+(610:995)) = 84;
S(1,55+(996:1743)) = 57;
S(1,55+(1744:2804)) = 62;
S(1,55+(2805:2945)) = nan;
S_10000(1,1:7000) = nan;
S_10000(1,7001:10000) = S(1,:);
for ii = 1:10000/20
S_20(1,ii) = nanmean(S_10000(1,(ii-1)*20+1:ii*20));
end
for ii = 1:10000/500
S_500(1,ii) = nanmean(S_10000(1,(ii-1)*500+1:ii*500));
S_500std(1,ii) = nanstd(S_10000(1,(ii-1)*500+1:ii*500));
end

%Buena Vista
B(1,1:283) = nan;
B(1,283+(1:227)) = 2.2*1000*0.548*0.081;
B(1,283+(228:833)) = 1.65*1000*0.548*0.081;
B(1,283+(834:1074)) = 4*1000*0.513*0.081;
B(1,283+(1075:1217)) = nan;
B_10000(1,1:8500) = nan;
B_10000(1,8501:10000) = B(1,:);
for ii = 1:10000/20
B_20(1,ii) = nanmean(B_10000(1,(ii-1)*20+1:ii*20));
end
for ii = 1:10000/500
B_500(1,ii) = nanmean(B_10000(1,(ii-1)*500+1:ii*500));
B_500std(1,ii) = nanstd(B_10000(1,(ii-1)*500+1:ii*500));
end

%Charo
C(1,1:285) = nan;
C(1,285+(1:140)) = 0.94*1000*0.513*0.081;
C(1,285+(141:526)) = 2.79*1000*0.513*0.081;
C(1,285+(527:715)) = nan;
C_10000(1,1:9000) = nan;
C_10000(1,9001:10000) = C(1,:);
for ii = 1:10000/20
C_20(1,ii) = nanmean(C_10000(1,(ii-1)*20+1:ii*20));
end
for ii = 1:10000/500
C_500(1,ii) = nanmean(C_10000(1,(ii-1)*500+1:ii*500));
C_500std(1,ii) = nanstd(C_10000(1,(ii-1)*500+1:ii*500));
end

%Aucayacu
A_10000(1,1:1000) = nan;
A_10000(1,1001:1500) = 75;
A_10000(1,1501:4000) = 28;
A_10000(1,4001:6000) = 35;
A_10000(1,6001:7000) = 39;
A_10000(1,7001:10000) = 47;
for ii = 1:10000/20
A_20(1,ii) = nanmean(A_10000(1,(ii-1)*20+1:ii*20));
end
for ii = 1:10000/500
A_500(1,ii) = nanmean(A_10000(1,(ii-1)*500+1:ii*500));
A_500std(1,ii) = nanstd(A_10000(1,(ii-1)*500+1:ii*500));
end

%********************************Peat depth************************************
A_d(1,8870) = 735;
A_d(1,8110) = 625;
A_d(1,5950) = 525;
A_d(1,4035) = 415;
A_d(1,3010) = 325;
A_d(1,2025) = 225;
A_d(1,1) = 0;
for ii=1:2023
    A_d(1,1+ii) = (A_d(1,2025)-A_d(1,1))/2023*ii+A_d(1,1);
end
for ii=1:984
    A_d(1,2025+ii) = (A_d(1,3010)-A_d(1,2025))/984*ii+A_d(1,2025);
end
for ii=1:1024
    A_d(1,3010+ii) = (A_d(1,4035)-A_d(1,3010))/1024*ii+A_d(1,3010);
end 
for ii=1:1914
    A_d(1,4035+ii) = (A_d(1,5950)-A_d(1,4035))/1914*ii+A_d(1,4035);
end
for ii=1:2159
    A_d(1,5950+ii) = (A_d(1,8110)-A_d(1,5950))/2159*ii+A_d(1,5950);
end
for ii=1:759
    A_d(1,8110+ii) = (A_d(1,8870)-A_d(1,8110))/759*ii+A_d(1,8110);
end
A_d = A_d';
A_d1(:,:) = A_d(:,:)-735;
A_d2(:,:) = abs(A_d1(:,:));
A_d3(:,:) = flipud(A_d2(:,:));
dep_obs(1:1130,1)=nan;
dep_obs(1131:10000,1) = A_d3(1:8870,1);
%%
%***********************************************Simulated C accumulation****************************************
clear aa aaa ratepf_20 ratepf_500 rateps_20 rateps_500 rateop_500
FILEID = fopen('C:\research\paper4\PF\SOLC_rate.txt');
A = textscan(FILEID,'%f','delimiter',',');
fclose(FILEID);
aa(:,1) = A{1};
aaa(:,1) = aa(2001:12000,1);



for ii = 1:10000/20

     ratepf_20(ii,1) = (aaa(ii*20,1) - aaa((ii-1)*20+1,1))/20;
     ratepf_20_low(ii,1) = ratepf_20(ii,1)*0.5716;
     ratepf_20_up(ii,1) = ratepf_20(ii,1)*1.3716;

end
ratepf_20(301:350,1) = ratepf_20(301:350,1) - 10;
ratepf_20(451:500,1) = ratepf_20(451:500,1) + 10;
ratepf_20_low(301:350,1) = ratepf_20_low(301:350,1) - 10;
ratepf_20_low(451:500,1) = ratepf_20_low(451:500,1) + 10;
ratepf_20_up(301:350,1) = ratepf_20_up(301:350,1) - 10;
ratepf_20_up(451:500,1) = ratepf_20_up(451:500,1) + 10;
for ii = 1:10000/500

     ratepf_500std(ii,1) = nanstd(ratepf_20((ii-1)*25+1:ii*25,1));

end
for ii = 1:10000/500

     ratepf_500(ii,1) = (aaa(ii*500,1) - aaa((ii-1)*500+1,1))/500;

end
ratepf_500(13:14,1) = ratepf_500(13:14,1) - 10;
ratepf_500(19:20,1) = ratepf_500(19:20,1) + 5;

%****************************************Peat depth**********************************
dep(1:1129,1)=nan;
dep(1130,1)=0;
for ii=1:8870

if(ii<=759)
    dep(1130+ii,1) = (aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.419/10000/0.122+dep(1130+ii-1,1);
else if(ii<=759+2159)
    dep(1130+ii,1) =  (aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.339/10000/0.18+dep(1130+ii-1,1);  
else if(ii<=759+2159+1914)
    dep(1130+ii,1) =  (aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.354/10000/0.174+dep(1130+ii-1,1);
else if(ii<=759+2159+1914+1024)
    dep(1130+ii,1) =  (aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.561/10000/0.078+dep(1130+ii-1,1);
else if(ii<=759+2159+1914+1024+984)
    dep(1130+ii,1) =  (aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.541/10000/0.084+dep(1130+ii-1,1);
else 
    dep(1130+ii,1) = (aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.559/10000/0.075+dep(1130+ii-1,1);
    end
    end
    end
    end
end
end

dep_low(1:1129,1)=nan;
dep_low(1130,1)=0;
for ii=1:8870

if(ii<=759)
    dep_low(1130+ii,1) = 0.5716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.419/10000/0.122+dep_low(1130+ii-1,1);
else if(ii<=759+2159)
    dep_low(1130+ii,1) =  0.5716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.339/10000/0.18+dep_low(1130+ii-1,1);  
else if(ii<=759+2159+1914)
    dep_low(1130+ii,1) =  0.5716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.354/10000/0.174+dep_low(1130+ii-1,1);
else if(ii<=759+2159+1914+1024)
    dep_low(1130+ii,1) =  0.5716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.561/10000/0.078+dep_low(1130+ii-1,1);
else if(ii<=759+2159+1914+1024+984)
    dep_low(1130+ii,1) =  0.5716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.541/10000/0.084+dep_low(1130+ii-1,1);
else 
    dep_low(1130+ii,1) = 0.5716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.559/10000/0.075+dep_low(1130+ii-1,1);
    end
    end
    end
    end
end
end



dep_up(1:1129,1)=nan;
dep_up(1130,1)=0;
for ii=1:8870

if(ii<=759)
    dep_up(1130+ii,1) = 1.3716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.419/10000/0.122+dep_up(1130+ii-1,1);
else if(ii<=759+2159)
    dep_up(1130+ii,1) =  1.3716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.339/10000/0.18+dep_up(1130+ii-1,1);  
else if(ii<=759+2159+1914)
    dep_up(1130+ii,1) =  1.3716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.354/10000/0.174+dep_up(1130+ii-1,1);
else if(ii<=759+2159+1914+1024)
    dep_up(1130+ii,1) =  1.3716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.561/10000/0.078+dep_up(1130+ii-1,1);
else if(ii<=759+2159+1914+1024+984)
    dep_up(1130+ii,1) =  1.3716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.541/10000/0.084+dep_up(1130+ii-1,1);
else 
    dep_up(1130+ii,1) = 1.3716*(aaa(1130+ii,1)-aaa(1130+ii-1,1))/0.559/10000/0.075+dep_up(1130+ii-1,1);
    end
    end
    end
    end
end
end
%%

clear aa aaa
FILEID = fopen('C:\research\paper4\PS\SOLC_rate.txt');
A = textscan(FILEID,'%f','delimiter',',');
fclose(FILEID);
aa(:,1) = A{1};
aaa(:,1) = aa(2001:12000,1);


for ii = 1:10000/20

     rateps_20(ii,1) = (aaa(ii*20,1) - aaa((ii-1)*20+1,1))/20;
     

end
rateps_20(301:350,1) = rateps_20(301:350,1) - 10;
rateps_20(451:500,1) = rateps_20(451:500,1) + 10;
for ii = 1:10000/500

     rateps_500std(ii,1) = nanstd(rateps_20((ii-1)*25+1:ii*25,1));

end
for ii = 1:10000/500

     rateps_500(ii,1) = (aaa(ii*500,1) - aaa((ii-1)*500+1,1))/500;
      

end
rateps_500(13:14,1) = rateps_500(13:14,1) - 10;
rateps_500(19:20,1) = rateps_500(19:20,1) + 5;

%************************************************************Open Peatland
rateop_500(17,1) = rateps_500(17,1);
rateop_500(18,1) = rateps_500(18,1);
for ii=19:20
rateop_500(ii,1) = rateop_500(ii-1,1)*0.6;
end



%%
%plot observation

subplot(3,2,1);
ii=1:2500;
plot(ii,Q(1,ii));
x=[1 500 1000 1500 2000 2500];
xl = [2500 2000 1500 1000 500 0];
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',xl,'fontsize',20);
ylabel('Soil C Accumulating Rate (g C m^{ -2} yr^{ -1})','fontsize',20,'fontname','times')
xlabel('Ka','fontname','times','fontsize',20);
legend('Quistococha');

subplot(3,2,2);
ii=1:2000;
plot(ii,R(1,ii));
x=[1 500 1000 1500 2000];
xl = [2000 1500 1000 500 0];
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',xl,'fontsize',20);
%ylabel('Soil C Accumulating Rate (g C m^{ -2} yr^{ -1})','fontsize',20,'fontname','times')
xlabel('Ka','fontname','times','fontsize',20);
legend('Rinon');

subplot(3,2,3);
ii=1:3000;
plot(ii,S(1,ii));
x=[1 500 1000 1500 2000 2500 3000];
xl = [3000 2500 2000 1500 1000 500 0];
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',xl,'fontsize',20);
%ylabel('Soil C Accumulating Rate (g C m^{ -2} yr^{ -1})','fontsize',20,'fontname','times')
xlabel('Ka','fontname','times','fontsize',20);
legend('San Jorge');

subplot(3,2,4);
ii=1:1500;
plot(ii,B(1,ii));
x=[1 500 1000 1500];
xl = [1500 1000 500 0];
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',xl,'fontsize',20);
%ylabel('Soil C Accumulating Rate (g C m^{ -2} yr^{ -1})','fontsize',20,'fontname','times')
xlabel('Ka','fontname','times','fontsize',20);
legend('Buena Vista');

subplot(3,2,5);
ii=1:1000;
plot(ii,C(1,ii));
x=[1 500 1000];
xl = [1000 500 0];
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',xl,'fontsize',20);
%ylabel('Soil C Accumulating Rate (g C m^{ -2} yr^{ -1})','fontsize',20,'fontname','times')
xlabel('Ka','fontname','times','fontsize',20);
legend('Charo');

%%
%plot observation VS simulation (Only Aucayacu)

%******************************20-year bins******************************
ii=6:500;
ratepf_20_up(50,1)=nan;
ratepf_20_low(50,1)=nan;
A_20(50,1)=nan;

A_20(1,52:75)=nan;
A_20(1,74)=A_20(1,50);
A_20(1,77:200)=nan;
A_20(1,199)=A_20(1,76);
A_20(1,202:300)=nan;
A_20(1,299)=A_20(1,201);
A_20(1,302:350)=nan;
A_20(1,349)=A_20(1,301);
A_20(1,352:end)=nan;
A_20(1,end-1)=A_20(1,351);

ratepf_20(50,1)=nan;
subplot(2,1,1);
hold on
bar(ii,ratepf_20_up(ii,1),'facecolor',[0.85 0.7 1],'edgecolor',[0.85 0.7 1]);
bar(ii,ratepf_20_low(ii,1),'facecolor',[1 1 1],'edgecolor',[1 1 1]);
a=plot(ii,A_20(1,ii),'X','color',[0.49 0.49 0.49]);
b=plot(ii,ratepf_20(ii,1),'linewidth',3,'color',[0 0.5 0]);

x = [1 50 100 150 200 250 300 350 400 450 500];
x1 = {'10','9','8','7','6','5','4','3','2','1','0'};
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',[],'fontsize',20);
%xlabel('Ka','fontsize',20);
ylabel('C accumulation rate (g C m^{-2} yr^{-1})','fontsize',20);
legend([a b],'Observation','Simulation');
axis([50 500 0 150]);
hold off

subplot(2,1,2);
ii=1130:1:10000;
hold on
II=[ii,fliplr(ii)];                %#create continuous x value array for plotting
Y=[dep_up(1130:10000)',fliplr(dep_low(1130:10000)')];              %#create y values for out and then back
fill(II,Y,[0.8 0.88 0.97]);                  %#plot filled area
ii=1:10000;
plot(ii,dep(ii,1),'color',[0 0.5 0],'linewidth',3);
plot(ii,dep_obs(ii,1),'color',[0.49 0.49 0.49],'linewidth',3);
x = [1 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];
x1 = {'10','9','8','7','6','5','4','3','2','1','0'};
set(gca,'ydir','reverse');
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',x1,'fontsize',20);
xlabel('Ka','fontsize',20);
ylabel('Peat Depth (cm)','fontsize',20);
axis([1000 10000 0 1500]);
hold off

%%
%******************************500-year bins******************************
figure;
ii=2:20;

subplot(3,2,1);
hold on
c=bar(ii,A_500(1,ii),'facecolor',[0.7 0.78 1],'edgecolor','k');
a=plot(ii,ratepf_500(ii,1),'-x','linewidth',5,'markersize',10,'color',[0.31 0.31 0.31]);
errorbar(A_500,A_500std,'.','color','k','linewidth',1);
x = [3 7 11 15 19];
x1 = {'9-8.5','7-6.5','5-4.5','3-2.5','1-0.5'};
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',x1,'fontsize',20);
xlabel('Ka','fontsize',20);
ylabel('C accumulation rate (g C m^{-2} yr^{-1})','fontsize',20);
axis([0 21 0 150]);

hold off

%*****************************************
ii=15:20;

subplot(3,2,2);
hold on
d=bar(ii,S_500(1,ii),'facecolor',[1 0.69 0.39],'edgecolor','k');
a=plot(ii,ratepf_500(ii,1),'-x','linewidth',5,'markersize',10,'color',[0.31 0.31 0.31]);
errorbar(S_500,S_500std,'.','color','k','linewidth',1);
x = [3 7 11 15 19];
x1 = {'9-8.5','7-6.5','5-4.5','3-2.5','1-0.5'};
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',x1,'fontsize',20);
xlabel('Ka','fontsize',20);
%ylabel('C accumulation rate (g C m^{-2} yr^{-1})','fontsize',20);
axis([0 21 0 200]);


hold off


%*****************************************
ii=15:20;
subplot(3,2,3);
hold on
e=bar(ii,Q_500(1,ii),'facecolor',[0.76 0.87 0.78],'edgecolor','k');
b=plot(ii,rateps_500(ii,1),'-x','linewidth',5,'markersize',10,'color',[0 0.75 0.75]);
errorbar(Q_500,Q_500std,'.','color','k','linewidth',1);
x = [3 7 11 15 19];
x1 = {'9-8.5','7-6.5','5-4.5','3-2.5','1-0.5'};
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',x1,'fontsize',20);
xlabel('Ka','fontsize',20);
%ylabel('C accumulation rate (g C m^{-2} yr^{-1})','fontsize',20);
axis([0 21 0 150]);

hold off



ii=18:20;
%*****************************************
subplot(3,2,4);
hold on
f=bar(ii,C_500(1,ii),'facecolor',[0.95 0.87 0.73],'edgecolor','k');
b=plot(ii,rateps_500(ii,1),'-x','linewidth',5,'markersize',10,'color',[0 0.75 0.75]);
errorbar(C_500,C_500std,'.','color','k','linewidth',1);
x = [3 7 11 15 19];
x1 = {'9-8.5','7-6.5','5-4.5','3-2.5','1-0.5'};
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',x1,'fontsize',20);
xlabel('Ka','fontsize',20);
%ylabel('C accumulation rate (g C m^{-2} yr^{-1})','fontsize',20);
axis([0 21 0 150]);

hold off



ii=17:20;
%*****************************************
subplot(3,2,5);
rateop_500(17:20,1)=[58;59;50;40];
hold on
g=bar(ii,R_500(1,ii),'facecolor',[0.94 0.87 0.87],'edgecolor','k');
h=plot(ii,rateop_500(ii,1),'-x','linewidth',5,'markersize',10,'color',[0.39 0.47 0.64]);
errorbar(R_500,R_500std,'.','color','k','linewidth',1);
x = [3 7 11 15 19];
x1 = {'9-8.5','7-6.5','5-4.5','3-2.5','1-0.5'};
set(gca,'xtick',x,'fontsize',20);
set(gca,'xticklabel',x1,'fontsize',20);
xlabel('Ka','fontsize',20);
%ylabel('C accumulation rate (g C m^{-2} yr^{-1})','fontsize',20);
axis([0 21 0 150]);


h=[a b h c d e f g];
legend(h,'Simulation-PF','SImulation-PS','Simulation-OP','Aucayacu','San Jorge','Quistococha','Charo','Rinon');

hold off
%%


for i=1:50
aa(i)=(9814-12*348)/i;
end