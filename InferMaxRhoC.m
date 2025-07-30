% clc
% clear all

load('MONOCLONAL_DATA.mat');
choosedata='R500';

format short
tinit=9;
tfinal=48;
tspan=linspace(tinit,tfinal,(tfinal-tinit)/3 + 1);

%choose concentration
chooseconc=1;

%scale the data and drop NaNs
[meanmatnorm,concvec]=ScaleData(choosedata);
[meanmatnorm,isnanmat, tspanvec, tspanmat]=DropNaNsFn(concvec,meanmatnorm,tspan);

data=meanmatnorm(chooseconc,1:tspanvec(chooseconc));
tspan=tspanmat(chooseconc,1:tspanvec(chooseconc));

figure
plot(tspan,data)

%initial scaled cell count
IC=data(1);

%constraints for inverse problem 
Aeq=[]; 
beq=[];
lb=0; 
ub=[];

%initial guess for rho
rho0=0.5;

gammas=1;
weights=ones(size(tspan));

gls_optpar=rho0; %initial optpar 
old_gls_optpar=rho0; %initial old optpar 
tol=1e-4; %tolerance: keeps weights above a certain size  
maxits=2000;  
minits=10; 
partol=0.1; %another tolerance, convergence criteria  
parchange=100; %initialize parameter change, start w smth large  
oldparchange=100; %initialize old parameter change  
ii=1; %initialize number of iterations 

while ii<maxits & parchange > partol & oldparchange > partol | ii< minits 
errtomin=@(par)ErrFnFindRhoC(par, data,tspan,IC,weights);
options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxIterations',550); 
weights;
[gls_optpar,~,~,~]=fmincon(errtomin, rho0, [], [], Aeq, beq, lb, ub, [], options); 
converge_flag=1;
[~,weights]=ForwardFnFindRhoC(gls_optpar,IC,tspan);
weights=weights';
weights(weights<tol)=0; 
weights(weights>tol)=weights(weights>tol).^(-2*gammas); 
inds=old_gls_optpar>1e-10; 
weights=full(weights); 
parchange =1/(2)*sum((abs(gls_optpar(inds)-old_gls_optpar(inds))./old_gls_optpar(inds))); 
ii = ii+1; 
old_gls_optpar=gls_optpar; 
end 
disp('GLS Estimation') 
gls_optpar

figure
[t,c]=ForwardFnFindRhoC(gls_optpar,IC, tspan);
plot(tspan,data,'*','LineWidth',2)
hold on
plot(tspan,c,'LineWidth',2)
legend('Data','Soln. w/ Max. rho','Location','NorthWest')
set(gca,"FontSize",20)
xlabel('Time')
ylabel('Most Resistant R. Pop.')

clr = hsv(length(concvec));

%same figure but fitted curve has color corresponding to concentration 
figure
[t,c]=ForwardFnFindRhoC(gls_optpar,IC, tspan);
plot(tspan,data,'k*','LineWidth',2)
hold on
plot(tspan,c,'LineWidth',2,'Color',clr(chooseconc,:))
legend('Fastest Growing Cells','Max. \rho Growth Term Fit','Location','NorthWest')
set(gca,"FontSize",20)
xlabel('Time (Hours)')
ylabel('Most Resistant R. Pop.')