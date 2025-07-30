function [tmat,gmat,ttreatmat,tgrowthmat,sprobsmat] = SyntheticCyclesFn(sgrid,sprobs,ttreatment,tgrowth,y0,treatcycles,growthcycles,noisesize,rho,k)

%inputs: initial distribution, time vectors for growth and treatment, 
%sgrid vector, y0 for first cycle, number of growth cycles, number of treatment cycles, noise size,
%rho, k, y0, 

% sgrid=linspace(0,1,11);
% sprobs=DistFn2('Normal',sgrid,0,1);
% ttreatment=linspace(0,10,11);
% tgrowth=linspace(0,5,6);
% y0=0.2;
% treatcycles=4;
% growthcycles=3; %growth cycles must be either same as treatcycles or one fewer
% noisesize=0;
% rho=0.3;
% k=0.45;

%sprobs matrix
sprobsmat=zeros(length(sgrid),treatcycles);
%sprobsmat(:,1)=sprobs;

%matrix of treatment cycles
tmat=zeros(length(ttreatment),treatcycles);

%matrix of growth cycles
timestep=ttreatment(2) - ttreatment(1);
gmat=zeros(length(tgrowth),treatcycles);
ttreatmat=zeros(length(ttreatment),treatcycles);
ttreatmat(:,1)=ttreatment;
tgrowthmat=zeros(length(tgrowth),growthcycles);
tgrowthmat(:,1)=tgrowth + ttreatment(end) + timestep;

%set up time vectors for all the treatment and growth cycles
% ttreatmat=zeros(length(ttreatment),treatcycles);
% tgrowthmat=zeros(length(tgrowth),growthcycles);
for i=2:treatcycles
    % ttreatmat(:,i)=ttreatmat(:,i-1) + ttreatmat(end,i-1)+ (i-1)*tgrowth(end) + 2*(i-1)*timestep;
    ttreatmat(:,i)=tgrowthmat(end,i-1)+ timestep + ttreatment;
    if i<=growthcycles
        % tgrowthmat(:,i)=tgrowthmat(:,i-1) + ttreatmat(end,i) + 2*(i-1)*timestep + timestep;
        tgrowthmat(:,i)= ttreatmat(end,i) + timestep + tgrowth;
    end
end

ttreatmat
tgrowthmat
%%
%run one cycle of treatment and one cycle of growth
for i=1:treatcycles
    sprobsmat(:,i)=sprobs;
    [t, cmat,weightedsol] = RK4FunctionC(sgrid, sprobs, rho, k, y0, ttreatment)
    weightedsoldata=weightedsol.*(1 + noisesize*randn(size(weightedsol)));
    tmat(:,i)=weightedsoldata;
    sprobs=(sprobs./(sum(sprobs))).*cmat(end,:); %new sprobs
    y0=weightedsoldata(end); %new IC
     %if i<=growthcycles

        [t, cmatgrow,weightedsolgrow] =RK4FunctionC(sgrid, sprobs, rho, 0, y0, tgrowth)
       %[t, cmatgrow,weightedsolgrow] =RK4FunctionC(sgrid, sprobs, rho, 0, y0, tgrowth )
        weightedsolgrowdata=weightedsolgrow.*(1 + noisesize*randn(size(weightedsolgrow)));
        gmat(:,i)=weightedsolgrowdata;
        y0=weightedsolgrow(end);

    %end
end
    
   
    % ttreatmat(:,i+1)=ttreatment + (i)*ttreatment(end) + (i)*tgrowth(end);
    % tgrowthmat(:,i+1)=tgrowth + (i+1)*ttreatment(end) + (i)*tgrowth(end);
%end

% %set up time vectors for all the treatment and growth cycles
% ttreatmat=zeros(length(ttreatment),treatcycles);
% tgrowthmat=zeros(length(tgrowth),growthcycles);
% for i=1:length(treatcycles)
%     ttreatmat(:,i)=ttreatment + (i-1)*ttreatment(end) + (i-1)*tgrowth(end);
%     if i<=growthcycles
%         tgrowthmat=ttreatmat(end,i) + tgrowth;
%     end
% end

% figure
% plot(ttreatmat(:,1), tmat(:,1),'*')
% hold on
% plot(ttreatmat(:,2), tmat(:,2),'*')
% hold on
% plot(ttreatmat(:,3), tmat(:,3),'*')

figure
for i=1:treatcycles
    hold on
    plot(ttreatmat(:,i),tmat(:,i),'k.','MarkerSize',12)
    %hold on
    if i<=growthcycles
        plot(tgrowthmat(:,i), gmat(:,i),'r.','MarkerSize',12 )
        hold on
    end
end
ylim([0 1])
legend('On Treatment','Off Treatment')
xlabel('Time')
ylabel('Aggregated Tumor Volume')
set(gca,"FontSize",20)
end