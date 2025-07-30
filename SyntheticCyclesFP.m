noisesize=0.0;
%dog=(sprobs./(sum(sprobs))).*cmat(end,:)
%% run forward problem first cycle for initial dist d1. save synthetic data.
% sgrid=linspace(0,1,101);
% sprobs=DistFn2('TwoPoints',sgrid,0,1);
sgrid=linspace(0,1,11);
sprobs=DistFn2('Normal',sgrid,0,1);
rho=0.3;
%k=0.45;
k=0.45;
y0=0.2;
tspan_data=linspace(0,10,11);
tspan_data_end=tspan_data(end);
time_step=tspan_data(2) - tspan_data(1);
[t, cmat,weightedsoldata] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan_data)

% for j=1:length(sprobs)
%     plot(tspan_data,cmat(:,j))
%     hold on
% end

figure
stem(sgrid,sprobs,'k*','LineWidth',2,'MarkerSize',8)


weightedsoldata=weightedsoldata.*(1 + noisesize*randn(size(weightedsoldata)));
% figure
% plot(tspan_data,weightedsoldata,'*')
%% compute d2 based on cmat from run 1 
%find proportions of total
totalpop=sum(cmat(end,:));
%sprobs2=cmat(end,:)./totalpop;
sprobs2=(sprobs./(sum(sprobs))).*cmat(end,:)
% figure
% plot(sgrid,sprobs2,'*')
% hold on
% plot(sgrid,sprobs,'d')


%% grow all pops from run 1 tfinal pops at rhomax for a while. save that
%synthetic data. 
%talk to Erica and see what I'm actually growing here. 
%For now, grow just the weighted sum at rhomax. 
%clearvars -except weightedsoldata sprobs sprobs2 tspan_data tspan_data_end noisesize time_step
clearvars tspan_grow
%sgridgrow=0 %this means total sensitivity, which is probably wrong 
%sprobsgrow=1
sgridgrow=sgrid;
sprobsgrow=sprobs2;
rho=0.3
y0=weightedsoldata(end)
%tspan_grow=linspace(0,1,11) + tspan_data(end)
tspan_grow=linspace(0,5,6) ;
tspan_grow_end=tspan_grow(end);
%[t, cmatgrow,weightedsolgrow] = RK4RhoGrowthFn(sgridgrow, sprobsgrow, rho, y0, tspan_grow)
[t, cmatgrow,weightedsolgrow] =RK4FunctionC(sgridgrow, sprobsgrow, rho, 0, y0, tspan_grow)
y02=weightedsolgrow(end);

weightedsolgrow=weightedsolgrow.*(1 + noisesize*randn(size(weightedsolgrow)));

%% cycle 2: run fp on d2
%sgrid=linspace(0,1,11);
%sprobs=sprobs2;
rho=0.3;
k=0.45;
%k=0.55;
y0=y02;
%tspan_data2=linspace(0,2,11) + tspan_grow(end);
tspan_data2=linspace(0,10,11) ;
tspan_data2_end=tspan_data2(end);
[t, cmat2,weightedsoldata2] = RK4FunctionC(sgrid, sprobs2, rho, k, y0, tspan_data2)

% figure
% for j=1:length(sprobs2)
%     plot(tspan_data2,cmat2(:,j))
%     hold on
% end
% hold on
% plot(tspan_data2,weightedsoldata2,'*','MarkerSize',12)

%compute s3
totalpop2=sum(cmat2(end,:));
%sprobs3=cmat2(end,:)./totalpop2;
sprobs3=(sprobs2./(sum(sprobs2))).*cmat2(end,:);

weightedsoldata2=weightedsoldata2.*(1 + noisesize*randn(size(weightedsoldata2)));


%% grow2
% sgridgrow2=0; %this means total sensitivity, which is probably wrong 
% sprobsgrow2=1;
sgridgrow2=sgrid;
sprobsgrow2=sprobs3;
rho=0.3;
y0=weightedsoldata2(end);
%tspan_grow=linspace(0,1,11) + tspan_data(end)
tspan_grow2=linspace(0,5,6) ;
tspan_grow2_end=tspan_grow2(end);
%[t, cmatgrow2,weightedsolgrow2] = RK4RhoGrowthFn(sgridgrow2, sprobsgrow2, rho, y0, tspan_grow2)
[t, cmatgrow2,weightedsolgrow2] = RK4FunctionC(sgridgrow2, sprobsgrow2, rho, 0, y0, tspan_grow2)
y03=weightedsolgrow2(end);

weightedsolgrow2=weightedsolgrow2.*(1 + noisesize*randn(size(weightedsolgrow2)));

%% %% cycle 3: run fp on d3
%sgrid=linspace(0,1,11);
%sprobs=sprobs2;
rho=0.3;
k=0.45;
%k=0.65;
y0=y03;
%tspan_data2=linspace(0,2,11) + tspan_grow(end);
tspan_data3=linspace(0,10,11) ;
tspan_data3_end=tspan_data3(end);
[t, cmat3,weightedsoldata3] = RK4FunctionC(sgrid, sprobs3, rho, k, y0, tspan_data3)


%compute s3
totalpop3=sum(cmat3(end,:));
%sprobs4=cmat3(end,:)./totalpop3;
sprobs4=(sprobs3./(sum(sprobs3))).*cmat3(end,:)

weightedsoldata3=weightedsoldata3.*(1 + noisesize*randn(size(weightedsoldata3)));

%% growth 3
% sgridgrow3=0; %this means total sensitivity, which is probably wrong 
% sprobsgrow3=1;
sgridgrow2=sgrid;
sprobsgrow3=sprobs4;
rho=0.3;
y0=weightedsoldata3(end);
%tspan_grow=linspace(0,1,11) + tspan_data(end)
tspan_grow3=linspace(0,5,6) ;
tspan_grow3_end=tspan_grow3(end);
%[t, cmatgrow3,weightedsolgrow3] = RK4RhoGrowthFn(sgridgrow3, sprobsgrow3, rho, y0, tspan_grow3)
[t, cmatgrow3,weightedsolgrow3] = RK4FunctionC(sgridgrow, sprobsgrow3, rho, 0, y0, tspan_grow3)
y04=weightedsolgrow3(end);

weightedsolgrow3=weightedsolgrow3.*(1 + noisesize*randn(size(weightedsolgrow3)));

%% compute s4
%sgrid=linspace(0,1,11);
%sprobs=sprobs2;
rho=0.3;
k=0.45;
%k=0.65;
y0=y04;
%tspan_data2=linspace(0,2,11) + tspan_grow(end);
tspan_data4=linspace(0,10,11) ;
tspan_data4_end=tspan_data4(end);
[t, cmat4,weightedsoldata4] = RK4FunctionC(sgrid, sprobs4, rho, k, y0, tspan_data4)


%compute s4
totalpop4=sum(cmat4(end,:));
%sprobs4=cmat3(end,:)./totalpop3;
sprobs5=(sprobs4./(sum(sprobs4))).*cmat4(end,:)

weightedsoldata4=weightedsoldata4.*(1 + noisesize*randn(size(weightedsoldata4)));
%% 
%plot all together
figure
plot(tspan_data,weightedsoldata,'k.','MarkerSize',12)
hold on
plot(tspan_grow + tspan_data_end + time_step, weightedsolgrow,'r.','MarkerSize',12)
hold on
plot(tspan_data2 + tspan_data_end + tspan_grow_end + 2*time_step,weightedsoldata2,'k.','MarkerSize',12)
hold on
plot(tspan_grow2 + tspan_data2_end + tspan_data_end + tspan_grow_end + 3*time_step, weightedsolgrow2,'r.','MarkerSize',12)
hold on
plot(tspan_data3 + tspan_grow2_end + tspan_data2_end + tspan_data_end + tspan_grow_end + 4*time_step, weightedsoldata3,'k.','MarkerSize',12)
hold on
plot(tspan_grow3 + tspan_data3_end + tspan_grow2_end + tspan_data2_end + tspan_data_end + tspan_grow_end + 5*time_step, weightedsolgrow3,'r.','MarkerSize',12)
hold on
plot(tspan_data4 + tspan_grow3_end + tspan_data3_end + tspan_grow2_end + tspan_data2_end + tspan_data_end + tspan_grow_end + 6*time_step, weightedsoldata4,'k.','MarkerSize',12)
set(gca,"FontSize",20)
ylim([0 1])
xlabel('time (unit lol)')
ylabel('tumor volume')
legend('on treatment','off treatment','Location','southeast')

% 
% %%
% figure
% stem(sgrid,sprobs./sum(sprobs),'d--','MarkerSize',6,'LineWidth',4)
% hold on
% stem(sgrid,sprobs2./sum(sprobs2),'*--','MarkerSize',6,'LineWidth',3)
% hold on
% stem(sgrid,sprobs3./sum(sprobs3),'o--','MarkerSize',6,'LineWidth',2)
% hold on
% stem(sgrid,sprobs4./sum(sprobs4),'x--','MarkerSize',6,'LineWidth',1)
% hold on
% stem(sgrid,sprobs5./sum(sprobs5),'s--','MarkerSize',6,'LineWidth',1)
% legend('p(s) cycle 1 (normalized)','p(s) cycle 2 (normalized)','p(s) cycle 3 (normalized','p(s) cycle 4 (normalized)','p(s) cycle 5 (normalized)')
% set(gca,"FontSize",20)


%%
mean1=sprobs./(sum(sprobs))*sgrid';
mean2=sprobs2./(sum(sprobs2))*sgrid';
mean3=sprobs3./(sum(sprobs3))*sgrid';
mean4=sprobs4./(sum(sprobs4))*sgrid';
mean5=sprobs5./(sum(sprobs5))*sgrid';

figure
plot(sgrid,sprobs./sum(sprobs),'b-','MarkerSize',6,'LineWidth',5)
hold on
plot(sgrid,sprobs2./sum(sprobs2),'r-','MarkerSize',6,'LineWidth',4)
hold on
plot(sgrid,sprobs3./sum(sprobs3),'m-','MarkerSize',6,'LineWidth',3)
hold on
plot(sgrid,sprobs4./sum(sprobs4),'k-','MarkerSize',6,'LineWidth',2)
hold on
plot(sgrid,sprobs5./sum(sprobs5),'g-','MarkerSize',6,'LineWidth',1)
hold on
plot(sgrid,ones(size(sgrid))*mean1,'bd','LineWidth',2,'MarkerSize',8)
hold on
plot(sgrid,ones(size(sgrid))*mean2,'r*','LineWidth',2,'MarkerSize',8)
hold on
plot(sgrid,ones(size(sgrid))*mean3,'mo','LineWidth',2,'MarkerSize',8)
hold on
plot(sgrid,ones(size(sgrid))*mean4,'kx','LineWidth',2,'MarkerSize',8)
hold on
plot(sgrid,ones(size(sgrid))*mean5,'gs','LineWidth',2,'MarkerSize',8)

legend('p(s) cycle 1 (normalized)','p(s) cycle 2 (normalized)','p(s) cycle 3 (normalized','p(s) cycle 4 (normalized','p(s) cycle 5 (normalized)','mean 1','mean 2','mean 3','mean 4','mean 5')
set(gca,"FontSize",14)
xlabel('s')
ylabel('p(s)')
% %% mean sensitivity for each cycle 
% mean1=sprobs./(sum(sprobs))*sgrid';
% mean2=sprobs2./(sum(sprobs2)).*sgrid';
% mean3=sprobs3./(sum(sprobs3)).*sgrid';
% mean4=sprobs4./(sum(sprobs4)).*sgrid';
% 
% plot()

% %%
% figure
% for j=1:length(sprobs)
%     plot(tspan_data,cmat(:,j),'k-')
%     hold on
% end
% for j=1:length(sprobs)
%     plot(tspan_data2,cmat(:,j),'r--')
%     hold on
% end
% for j=1:length(sprobs)
%     plot(tspan_data3,cmat(:,j),'b*')
%     hold on
% end