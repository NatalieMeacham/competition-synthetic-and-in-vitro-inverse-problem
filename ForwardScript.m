%Solve forward problem for a given initial distribution, growth rate, death
%rate, and initial condition
%Generate figure 1 from paper 

%establish sensitivity distribution
a=0;
b=1;
points=11;
sgrid=linspace(a,b,points);
[sprobs] = DistFn2('Bigaussian',sgrid,a,b);

%other parameters for forward solution 
rho=1;
k=1.5;
y0=.2; 
tfinal=10;
tpoints=101;
tspan=linspace(0,tfinal,tpoints); 

%compute forward solution
%[t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan);
[t, cmat,weightedsol] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan);

%plot distribution
figure
for i=1:length(sgrid)
    hold on
    stem(sgrid(i), sprobs(i),'--o','LineWidth',3,'MarkerSize',12)
end
set(gca,"FontSize",20)
xlabel('Sensitivity to Treatment{\it s}')
ylabel('Proportion of Population')

%plot all c(t,s) curves
figure
for i=1:length(sgrid)
    hold on
    plot(t, cmat(:,i),'-','LineWidth',2)
    Legend{i}=strcat('{\it s=}',' ', num2str(sgrid(i)));
end
set(gca,"FontSize",20)
xlabel('Time (Hours)')
ylabel('Subpopulation Growth')
legend(Legend)
ylim([0,1])

%plot aggregated tumor volume
figure
plot(t,weightedsol,'k','LineWidth',2)
xlabel('Time (Hours)')
ylabel('Aggregated Volume')
ylim([0,1])
set(gca,"FontSize",20)