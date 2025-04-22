clear all
points=10;
disttype='Bigaussian';
rho=0.3;
k=0.45;
y0=0.2;
tfinal=10;
tpoints=25;
noisesize=0.02;
rpoints=10;
a=0;
b=1;
tspan=linspace(0,tfinal,tpoints);
sgrid=linspace(0,1,points);
sprobs = DistFn2(disttype,sgrid,a,b);


%create propdata
[t, cmat,weightedsol] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan); %%regular
%[t, cmat,weightedsol] = RK4FunctionC(sgrid, ones(size(sprobs)), rho, k,y0, tspan); %experiment 
[t, cmatR,~] = RK4FunctionC(linspace(0,1,rpoints), ones(1,rpoints), rho, k, y0, tspan);
propdata=weightedsol.*(1 + noisesize*randn(size(weightedsol)));

[gls_optpar,converge_flag,AIC_GLS,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFnC(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpoints,propdata)

figure
plot(sgrid,sprobs)
hold on
plot(rsgrid,gls_optpar)

figure
plot(t,propdata,'*')
hold on
plot(t,cmatR*gls_optpar)
legend('propdata','recovered fit')

figure
for i=1:points
    plot(t,cmat(:,i),'b')
    hold on
end
for j=1:rpoints
    plot(t,cmatR(:,j),'r--')
    hold on
end