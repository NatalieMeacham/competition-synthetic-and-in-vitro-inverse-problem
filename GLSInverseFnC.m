function[gls_optpar,converge_flag,AIC_GLS,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFnC(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpoints,propdata)

%FORWARD PROBLEM
%Create synthetic data
 a=0; 
 b=1;
 sgrid=linspace(a,b,points);
% 
% %Create original distribution
[sprobs] = DistFn2(disttype,sgrid,a,b);

%Compute tspan for ode45 to use
tspan=linspace(0,tfinal,tpoints); 

%Other parameters
% rho=1;
% k=1.5;
% y0=.2; 

%Forward function to get weightedsol
[t, cmat,weightedsol] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan);

%Create synthetic data by adding some noise to weightedsol:
%noisesize=0.1;
%propdata = weightedsol.*(1+noisesize*rand(size(weightedsol)));

%PARAMETER ESTIMATION
rsgrid = linspace(a,b,rpoints);

%Set up constraints for constrained optimization (ie s must sum to 1 and have nonneg vals)
%Note that these are constraints for what we're trying to recover
%Must sum to 1:
Aeq=zeros(1,rpoints); %this chunk does not want to be fn inputs 
Aeq(1,:)=1;
beq=ones(1,1);
%Only nonneg vals:
lb=zeros(rpoints,1); 
ub=zeros(rpoints,1);
ub(1:rpoints)=1;

%Define initial s
%s0=rand(rpoints,1); %random s0

%uniform s0
s0=ones(rpoints,1); %to get a uniform dist.
s0=s0/sum(s0); %normalize
init_guess=s0';

%Weights for use in PEerrorfn 
gammas=1;
weights=ones(size(weightedsol)); %doesn't need to be input

gls_optpar=init_guess; %initial optpar 
old_gls_optpar=init_guess; %initial old optpar 
tol=1e-4; %tolerance: keeps weights above a certain size  
maxits=2000;  
minits=10; 
partol=0.1; %another tolerance, convergence criteria  
parchange=100; %initialize parameter change, start w smth large  
oldparchange=100; %initialize old parameter change  
ii=1; %initialize number of iterations  

%while ii<maxits && parchange > partol && oldparchange > partol || ii< minits 
while ii<maxits & parchange > partol & oldparchange > partol | ii< minits 
%gls_error_estimate = @(par)gls_formulation(par,prop_data,ts,y0,weights); 
%[t, gencmatC,~] = RK4FunctionC(rsgrid, ones(size(rsgrid)), rho, k, y0, tspan); %former 
[t, gencmatC,~] = RK4FunctionC(rsgrid, s0', rho, k, y0, tspan); %experimetal
errtomin=@(par)PEerrorfn(par,gencmatC,propdata,weights);
options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxIterations',550); 
%gls_optpar=FMINCON (set up Aeq etc.) 
[gls_optpar,~,converge_flag,~]=fmincon(errtomin, s0, [], [], Aeq, beq, lb, ub, [], options);
%pause
[~, ~,weights] = RK4FunctionC(rsgrid, gls_optpar', rho, k, y0, tspan); 
weights(weights<tol)=0; 
weights(weights>tol)=weights(weights>tol).^(-2*gammas); 
inds=old_gls_optpar>1e-10; 
weights=full(weights); 
parchange =1/(2)*sum((abs(gls_optpar(inds)-old_gls_optpar(inds))./old_gls_optpar(inds))); 
ii = ii+1; 
old_gls_optpar=gls_optpar; 
s0=gls_optpar;
end 
disp('GLS Estimation') 
gls_optpar 
 
%solve w optpar 
%compute residuals between data and gls_sol 

%[t,gencmat,~] = ForwardFunction(rsgrid, ones(size(rsgrid)), rho, k, y0, tspan); %I
%don't actually need a general matrix do i? aside from the name in error
%fn? Yes I think so bc it has different # of grid points
%[t, gencmatC,~] = ForwardFunctionN(rsgrid, ones(size(rsgrid)), rho, k, y0, tspan);
%note here I replaced sprobs with ones(size(rsgrid)) which I think is
%correct but not completely sure

%Create error function to minimize with fmincon
%errtomin=@(s)PEerrorfn(s,gencmatC,constdata,weights); %define error function where r is unknown aka anonymous fn over r

%Options for fmincon
%options = optimoptions('fmincon') %default options for fmincon
%options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxIterations',550); %matlab gives 1500

%Minimize errtomin to obtain optweight, convergence, etc.
%[optweight,~,converge_flag,~]=fmincon(errtomin, s0, [], [], Aeq, beq, lb, ub, [], options);

%Save optweight to a file in matlab drive - check where files are landing
%bc might not be in spring file
%savename = strcat('RS','rp',num2str(rpoints),disttype);
%save(savename,'optweight')

%AIC SECTION
%take error where s is given by optweight for AIC
%[current_err,~] = errorfunc_discrete_sparse(optweight,fullsol,agg_sol,weights);
[finalerr,~]=PEerrorfn(gls_optpar, gencmatC, propdata,weights)
gls_optpar
weights


%Compute AIC score for the given points and rpoints here:
%AIC_OLS = (points)*log(finalerr/points) + (points)*log((2*pi)+1) + 2*(rpoints + 1);
%AIC_OLS = (tpoints)*log(finalerr/tpoints) + 2*(rpoints + 1);
AIC_GLS = (tpoints)*log(sum((weights).*((propdata-(gencmatC*gls_optpar)').^2))/tpoints) + 2*(rpoints + 1);
%AIC_GLS = (tpoints)*log(sum((weights).*((propdata-(gencmatC*gls_optpar)').^2))/tpoints) + 2*(rpoints + 1);
%AIC_GLS = (tpoints)*log(sum((weights.^(-2)).*((propdata-(gencmatC*gls_optpar)').^2))/tpoints) + 2*(rpoints + 1) + 2*(rpoints+1)*(rpoints+2)/(tpoints - rpoints);

%NOTE: Banks hints that regular AIC should be used only if the sample size
%is at least 40x the number of estimated parameters. 

% %compare output of fmincon with original dist
% figure
% yyaxis left
% plot(rsgrid,optweight,'LineWidth',2)
% hold on
% yyaxis right
% plot(sgrid,sprobs,'LineWidth',2)
% legend('recovered dist.','original dist.')
% xlabel('sgrid')
% ylabel('s distribution')
% %figurename = strcat('F','rp',num2str(rpoints),'sp',num2str(points));
% %saveas(gcf,figurename)

%%compare weighted csol from optweight with original weightedsol
%[~,~,sweightedsol] = ForwardFunction(rsgrid, optweight', rho, k, y0, tspan);
%figure 
%plot(t,sweightedsol,'LineWidth',2)
%hold on
%plot(t,weightedsol','*')
%hold off
%legend('Soln. w. optweight','Soln. w. sprobs')
%%keyboard
end