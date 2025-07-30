function[err]=ErrFnFindRhoC(rho, data,tspan,IC,weights)

%solve forward problem
%[~,c]=ode45(@(t,c) rho*c*(1-c),tspan,IC );
[~,~,c] = ForwardFnFindRhoC(rho, IC, tspan);

%compute error
err=sum((weights).*((c - data)).^2);