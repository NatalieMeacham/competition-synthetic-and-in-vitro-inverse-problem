function [AICminvec, trymatoptrsgridS, optmatS, cdfRmatS, wsolmat,evvec,concvecS]=GLSAllDataInclNaNsC(choosedata,figures)

%Inverse problem built around the needs of the in vitro datasets
%Note that you can choose to return figures or not as 'y' or 'n'

%INPUTS
    %choosedata: choice of dataset
    %figures: yes or no for whether this function returns figures

%OUTPUTS
    %AICminvec: minimum AIC value for each concentration
    %trymatoptrsgridS: matrix of grids for the recovered distribution for each dosage
    %optmatS: matrix of optimal distribution for each dosage
    %cdfRmatS: matrix of recovered CDFs by dosage
    %wsolmat: aggregated tumor volume from each dosage's optimal
    %recovered distribution
    %evvec: mean recovered sensitivity for each dosage
    %concvecS: vector of relevant concentration numbers 

format long

%rho and k values for mixture datasets
kmax=0.0077; 
kmin=0;
% rhomax=0.0692; %old version
rhomax=0.0637;
%rhomin=0.0375; %old version
rhomin=0.0233;
%rho and k values for monoclonal datasets
if choosedata=='R250'
    %rhomin=0.0241;
    rhomin=0.0233;
    %rhomax=0.0692; %old model
    rhomax=0.0637;
    kmin=0;
    kmax=0.0077; 
    %kmax=0.0637;
elseif choosedata=='R500'
    %rhomin=0.0375; %old model
    rhomin=0.0365;
    %rhomax=0.0682;  %old model
    rhomax=0.0634;
    kmin=0;
    kmax=0.0077; 
    %kmax=0.0634;
elseif choosedata=='S500'
    %rhomin=0.0241; %old version
    rhomin=0.0233;
    %rhomax=0.0692; %old version
    rhomax=0.0637;
    kmin=0;
    kmax=0.0019; 
    %kmax=0.0637;
elseif choosedata=='S100'
    %rhomin=0.0241;  %old version
    rhomin=0.0233;
    %rhomax=0.0692; %old version
    rhomax=0.0637;
    kmin=0;
    kmax=0.0077; 
    %kmax=0.0637;
else
    disp('??')
end
%kmax=rhomax;

dosevec= [0, 0.0260,0.0536, 0.1094, 0.2244, 0.3410,0.4603,1.1908,2.4430,3.7133,5.0119];
dosevecshort= [0.00, 0.03,0.05, 0.12, 0.22, 0.34,0.46,1.19,2.44,3.71,5.01];
dosevecstring= ["0.00", "0.03","0.05", "0.12", "0.22", "0.34","0.46","1.19","2.44","3.71","5.01"];

%scale k by dosage
kvec=kmax.*(dosevec./dosevec(11)); 

%create vector of max rho value
rhovec=rhomax*ones(length(dosevec)); 

tinit=9;
tfinal=48;
tspan=linspace(tinit,tfinal,(tfinal-tinit)/3 + 1);

[meanmatnorm,concvec,errmat]=ScaleData(choosedata);

[meanmatnorm,isnanmat, tspanvec, tspanmat,concvecA]=DropNaNsFn(concvec,meanmatnorm,tspan);

%find AIC vals for each concentration for each number of rpoints
rpointsvec=4:1:30;
AICmat=zeros(length(concvecA),length(rpointsvec));
for a=1:length(concvecA)
    data=meanmatnorm(a,1:tspanvec(a));
    y0=data(1);
    rho=rhovec(concvecA(a));
    k=kvec(concvecA(a));
    tspanS=tspanmat(a,1:tspanvec(a));
    for i=1:length(rpointsvec)
    [~,~,AIC_GLS,~,~]=GLSInverseFnCData(rho,k,y0,tspanS,rpointsvec(i),data);
    AICmat(a,i)=AIC_GLS;
    end
end

%create a matrix that fills in the optweight corresponding to the lowest
%AIC value for each concentration
concnum=length(concvecA);
AICminvec=zeros(concnum,1);
AICcheck=zeros(concnum,1);
optmat=zeros(concnum,length(rpointsvec));
for a=1:length(concvecA)
    data=meanmatnorm(a,1:tspanvec(a));
    tspanS=tspanmat(a,1:tspanvec(a));
    rho=rhovec(concvecA(a));
    k=kvec(concvecA(a));
    [M,I]=min(AICmat(a,:));
    AICminvec(a)=rpointsvec(I);
    y0=data(1);
    [optweight,~,~,~,~]=GLSInverseFnCData(rho,k,y0,tspanS,AICminvec(a),data);
    optmat(a,1:AICminvec(a))=optweight;
end

concvecS=concvecA;

%create legend for plots
for a=1:length(concvecS)
    for iter=1:length(concvecS)
        Legend{iter}=strcat('Conc.',' ', num2str(concvecS(iter)));
        Legend2{iter}=strcat('d=',num2str(dosevecshort(concvecS(iter))),'  \muM');
        Legend3{iter}=strcat('{\it d}=',dosevecstring(concvecS(iter)));
    end
end

clr = hsv(length(concvecS));

%create matrices of optimal grids and optimal dists 
trymatoptrsgrid=zeros(length(concvecS),max(AICminvec));
for a=1:length(concvecS)
    optrpoints=AICminvec(a);
    optrsgrid=linspace(0,1,optrpoints);
    trymatoptrsgrid(a,1:optrpoints)=optrsgrid;
end

%create matrix and compute cdfs for all dists
cdfRmat=zeros(length(concvecS),max(AICminvec));
for a=1:length(concvecS)
    optrpoints=AICminvec(a);
    optrsgrid=linspace(0,1,optrpoints);
    partialsumR=0;
    for i=1:optrpoints
        partialsumR=sum(optmat(a,1:i));
        cdfRmat(a,i)=partialsumR;
    end
end

trymatoptrsgridS=trymatoptrsgrid;
optmatS=optmat; 
AICminvecS=AICminvec;
cdfRmatS=cdfRmat; 
meanmatnormS=meanmatnorm;

wsolmat=zeros(length(concvecS),length(tspan));
for a=1:length(concvecS)
    rho=rhovec(concvecA(a));
    k=kvec(concvecA(a));
    data=meanmatnorm(a,1:tspanvec(a));
    tspanS=tspanmat(a,1:tspanvec(a));
    y0=data(1);
    optrpoints=AICminvecS(a);
    optrsgrid=linspace(0,1,optrpoints);
    optweight=optmatS(a,1:AICminvecS(a));
    [~, ~,weightedsol] = RK4FnCData(optrsgrid, optweight, rho, k, y0, tspanS);
    wsolmat(a,1:length(tspanS))=weightedsol;
end

if figures == 'y'

%plot PDF/PMFs
figure
for a=1:length(concvecS)
    hold on
    stem(trymatoptrsgridS(a,1:AICminvecS(a)),optmatS(a,1:AICminvecS(a)),'--o','Color',clr(a,:),'MarkerFaceColor',clr(a,:),'MarkerSize',8,'LineWidth',2)    
end

xlabel('Sensitivity to Treatment {\it s}')
ylabel('Recovered Proportion of Population')
legend(Legend, 'FontSize',12)
legend(Legend2, 'FontSize',12)
legend(Legend3, 'FontSize',12,'Location','northwest')
set(gca,"FontSize",20)
saveas(gcf,strcat('Figures/',choosedata,'PDF.fig'));
saveas(gcf,strcat('Figures/',choosedata,'PDF.png'));

%plot CDFs
figure
for a=1:length(concvecS)
    hold on
    clr = hsv(length(concvecS));
    X=trymatoptrsgridS(a,1:AICminvecS(a));
    Y=cdfRmatS(a,1:AICminvecS(a));
    stairs(X,Y,'--o','Color',clr(a,:),'MarkerFaceColor',clr(a,:),'MarkerSize',8,'LineWidth',2)
    ylim([0 1])
    xlabel('Sensitivity to Treatment {\it s}')
ylabel('Cumulative Recovered Proportion')
legend(Legend3, 'FontSize',12,'Location','northwest')
set(gca,"FontSize",20)
end
saveas(gcf,strcat('Figures/',choosedata,'CDF.fig'));
saveas(gcf,strcat('Figures/',choosedata,'CDF.png'));

%plot fits
figure
for a=1:length(concvecS)
    rho=rhovec(concvecA(a));
     k=kvec(concvecA(a));
     data=meanmatnorm(a,1:tspanvec(a));
    tspanS=tspanmat(a,1:tspanvec(a));
     y0=data(1);
     optrpoints=AICminvecS(a);
     optrsgrid=linspace(0,1,optrpoints);
     optweight=optmatS(a,1:AICminvecS(a));
   [~, ~,weightedsol] = RK4FnCData(optrsgrid, optweight, rho, k, y0, tspanS);
    clr = hsv(length(concvecS));    
    errorbar(tspanS,data,errmat(a,1:tspanvec(a)),'.','MarkerSize',8,'LineWidth',2,'Color',clr(a,:),'HandleVisibility','off')
    hold on
    plot(tspanS,weightedsol,'LineWidth',2,'Color',clr(a,:))
    hold on
end
xlabel('Time (Hours)')
    ylabel('Tumor Growth Data and Model Fit')
    legend(Legend3, 'FontSize',12,'Location','northwest')
    set(gca,"FontSize",20)
saveas(gcf,strcat('Figures/',choosedata,'Fit.fig'));
saveas(gcf,strcat('Figures/',choosedata,'Fit.png'));
else
    disp('no figures')
end

%compute expected value of the recovered distribution
evvec=zeros(length(concvecS),1);
for a=1:length(concvecS)
    optrpoints=AICminvecS(a);
    optrsgrid=linspace(0,1,optrpoints);
    optweight=optmatS(a,1:AICminvecS(a));
    evvec(a)=sum(optrsgrid.*optweight);
end