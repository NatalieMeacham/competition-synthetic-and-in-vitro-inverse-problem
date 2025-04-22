a=0;
b=1;
points=11;
sgrid=linspace(a,b,points);

sprobs=DistFn2('Uniform',sgrid,a,b);

rho=0.3;
k=0.45;
y0=0.2;
tfinal=10;
tpoints=25;
tspan=linspace(0,tfinal,tpoints);

[t, cmat,weightedsol,uvec] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan)

figure
plot(t,uvec,'*')
hold on
plot(t,weightedsol)
legend('uvec','weightedsol')