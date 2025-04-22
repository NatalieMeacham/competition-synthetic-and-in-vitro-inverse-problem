function [t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan)
    %Solve ODE for each point on sgrid:
    
    cmat=zeros(length(tspan),length(sgrid));
    %figure
    for i=1:length(sgrid)
        [t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) -k*sgrid(i)*csol, tspan, y0);
        %[t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) - k*sgrid(i)*csol*(1 - csol), tspan, y0); %logistic version
        %[t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) - k*rho*sgrid(i)*csol, tspan, y0); %version where death term is scaled by growth rate
        %fill in matrix for csol values over time and different s values
        cmat(:,i)=csol;
    end

    %weighted sum with normalized probabilities
    weightedsol=(sprobs./(sum(sprobs)))*cmat';
 end 