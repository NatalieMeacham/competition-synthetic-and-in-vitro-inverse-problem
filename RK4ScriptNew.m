a=0;
b=1;
points=11; %Note to self: the pointval*100 is the problem for different numbers of points
sgrid=linspace(a,b,points);

%[sprobs] = DistFn(['OnePoint'],sgrid,a,b);
[sprobs] = DistFn2('Normal',sgrid,a,b);

%Other parameters for forward solution 
rho=.3;
k=.45;
y0=.2; %note c is normalized so needs to start between 0 and 1
tfinal=50;
tpoints=26;
tspan=linspace(0,tfinal,tpoints); %this is so we can initialize matrix with known number of rows 

%[t, cmat] = ForwardFunction(sgrid, sprobs, rho, k, y0, tspan)
[t, cmat,weightedsol] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan);

%%
%run forward solver for time spans with 100, 25, 10, and 5 time steps 
a=0;
b=1;
points=11; %Note to self: the pointval*100 is the problem for different numbers of points
sgrid=linspace(a,b,points);

%[sprobs] = DistFn(['OnePoint'],sgrid,a,b);
[sprobs] = DistFn2('Normal',sgrid,a,b);

%Other parameters for forward solution 
rho=.3
k=.45;
y0=.2; %note c is normalized so needs to start between 0 and 1

t1=101;
t2=26;
t3=11;
t4=6;

tfinal=10;
tspan1=linspace(0,tfinal,t1)
tspan2=linspace(0,tfinal,t2)
tspan3=linspace(0,tfinal,t3)
tspan4=linspace(0,tfinal,t4)

[tvec1, cmat1,weightedsol1] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan1);
[tvec2, cmat2,weightedsol2] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan2);
[tvec3, cmat3,weightedsol3] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan3);
[tvec4, cmat4,weightedsol4] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan4);

figure
plot(tvec1,weightedsol1)
hold on
plot(tvec2,weightedsol2)
hold on
plot(tvec3,weightedsol3)
hold on
plot(tvec4,weightedsol4)
legend('1','2','3','4')
ylim([0 1])
%%
figure
for i=1:points
    plot(t,cmat(:,i))
    hold on
end
hold on
plot(t,weightedsol,'*')
%%
    figure
    for i=1:length(sgrid)
        hold on
        stem(sgrid(i), sprobs(i),'--o','LineWidth',3,'MarkerSize',12)
    end
    set(gca,"FontSize",20)
    %title('Proportions of Different s Values')
    xlabel('Sensitivity to Treatment (s)')
    ylabel('Proportion of Population')
    %ylim([0,1])
%%
% %plot one curve
% figure
% plot(t,cmat(1,:))

    % %plot all curves
    % figure
    % for i=1:length(sgrid)
    %     hold on
    %     plot(t, cmat(:,i),'-','LineWidth',2)
    %     Legend{i}=strcat('s=',' ', num2str(sgrid(i)));
    % end
    % set(gca,"FontSize",20)
    % %title('Proportions of Different s Values')
    % xlabel('Time')
    % ylabel('Subpopulation Growth')
    % legend(Legend)
    % ylim([0,1])

    %plot weightedsol
    figure
    plot(t,weightedsol,'k','LineWidth',2)
    xlabel('Time')
    ylabel('Aggregated Volume')
    ylim([0,1])
    set(gca,"FontSize",20)

   
