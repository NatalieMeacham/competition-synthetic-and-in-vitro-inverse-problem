function timestepplots(sprobs100,optweightfromAIC100,optweightfromAIC25,optweightfromAIC10,optweightfromAIC5,weightedsol100,weightedsol25,weightedsol10,weightedsol5,Sweightedsol100,Sweightedsol25,Sweightedsol10,Sweightedsol5,disttype,t,weightedsol,noisesize)
points = 101;
noisestr=string(noisesize);

%calculate the cdfs
if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1 
    cdfS100 = cumsum(sprobs100);
elseif strcmp(disttype,'Normal') == 1 || strcmp(disttype,'Uniform') ==1 || strcmp(disttype,'Bigaussian') ==1 
    cdfS100=cumtrapz(linspace(0,1,length(sprobs100)),sprobs100);
end

%plot cdfs
figure
stairs(linspace(0,1,length(optweightfromAIC100)),cumsum(optweightfromAIC100),'-','MarkerSize',8,'LineWidth',2,'Color','#0072BD')
hold on
stairs(linspace(0,1,length(optweightfromAIC25)),cumsum(optweightfromAIC25),'-','MarkerSize',8,'LineWidth',2,'Color','#6B9C28')
hold on
stairs(linspace(0,1,length(optweightfromAIC10)),cumsum(optweightfromAIC10),'-','MarkerSize',8,'LineWidth',2,'Color','#E6AB1A')
hold on
stairs(linspace(0,1,length(optweightfromAIC5)),cumsum(optweightfromAIC5),'-','MarkerSize',8,'LineWidth',2,'Color','k')
    ylim([0 1])
    hold on
    plot(linspace(0,1,points),cdfS100,'-','MarkerSize',6,'LineWidth',2,'Color','#A2142F')
    xlabel('Sensitivity to Treatment {\it s}')
%ylabel('Cumulative Recovered Prop. of Pop.') %original version
ylabel('Cumulative Recovered Proportion')
%legend(Legend, 'FontSize',12,'Location','northwest')
legend('Recovered, N_t=100','Recovered, N_t=25','Recovered, N_t=10','Recovered, N_t=5','Original','Location','best','FontSize',12)
set(gca,"FontSize",20)
CDFfiglabel=strcat(disttype,'N',noisestr,'_timepoints_CDF','.jpg');
saveas(gcf,CDFfiglabel);
CDFfiglabel=strcat(disttype,'N',noisestr,'_timepoints_CDF','.fig');
saveas(gcf,CDFfiglabel);


%plot fits
figure
plot(linspace(0,50,101),weightedsol100,'d','LineWidth', 2,'MarkerSize',8,'Color','#0072BD')
hold on
plot(linspace(0,50,101),Sweightedsol100,'LineWidth',2,'Color','#0072BD')
hold on
plot(linspace(0,50,26),weightedsol25,'d','LineWidth', 2,'MarkerSize',8,'Color','#6B9C28')
hold on
plot(linspace(0,50,26),Sweightedsol25,'LineWidth',2,'Color','#6B9C28')
hold on
plot(linspace(0,50,11),weightedsol10,'d','LineWidth', 2,'MarkerSize',8,'Color','#E6AB1A')
hold on
plot(linspace(0,50,11),Sweightedsol10,'LineWidth',2,'Color','#E6AB1A')
hold on
plot(linspace(0,50,6),weightedsol5,'d','LineWidth', 2,'MarkerSize',8,'Color','k')
hold on
plot(linspace(0,50,6),Sweightedsol5,'LineWidth',2,'Color','k')
hold on
plot(t',weightedsol,'-','LineWidth',2,'Color','#A2142F')
legend('Synthetic Data, N_t=100','Recovered Curve, N_t=100','Synthetic Data, N_t=25','Recovered Curve, N_t=25','Synthetic Data, N_t=10','Recovered Curve, N_t=10','Synthetic Data, N_t=5','Recovered Curve, N_t=5','Original Curve','Location','best','FontSize',12)
set(gca,"FontSize",20)
ylabel('Aggregated Tumor Volume')
xlabel('Time')
ylim([0 1])
Fitfiglabel=strcat(disttype,'N',noisestr,'_timepoints_Fit','.jpg');
saveas(gcf,Fitfiglabel);
Fitfiglabel=strcat(disttype,'N',noisestr,'_timepoints_Fit','.fig');
saveas(gcf,Fitfiglabel);



r100to25=length(optweightfromAIC25)/length(optweightfromAIC100);
r100to10=length(optweightfromAIC10)/length(optweightfromAIC100);
r100to5=length(optweightfromAIC5)/length(optweightfromAIC100);

%plot pmfs
if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1 
    figure
  %  yyaxis left
    stem(linspace(0,1,length(optweightfromAIC100)),optweightfromAIC100,'--square','LineWidth',2,'MarkerSize',10,'Color','#0072BD')
    hold on
    stem(linspace(0,1,length(optweightfromAIC25)),optweightfromAIC25,'--*','LineWidth',2,'MarkerSize',10,'Color','#6B9C28')
    hold on
    stem(linspace(0,1,length(optweightfromAIC10)),optweightfromAIC10,'--v','LineWidth',2,'MarkerSize',10,'Color','#E6AB1A')
    hold on
    stem(linspace(0,1,length(optweightfromAIC5)),optweightfromAIC5,'--x','LineWidth',2,'MarkerSize',10,'Color','k')
    hold on
    ylabel('Recovered Proportion of Population')
    %yyaxis right
    stem(linspace(0,1,101),sprobs100,'--o','LineWidth',2,'MarkerSize',10,'Color','#A2142F')
    %ylabel('Original Proportion of Population')
    ylimval=max([max(optweightfromAIC100),max(optweightfromAIC25),max(optweightfromAIC10),max(optweightfromAIC5),max(sprobs100)]);
    ylim([0 ylimval])
    legend('Recovered, N_t=100','Recovered, N_t=25','Recovered, N_t=10','Recovered, N_t=5','Original','Location','best','FontSize',12)
    set(gca,"FontSize",20)
    xlabel('Sensitivity to Treatment {\it s}')
elseif strcmp(disttype,'Normal') == 1 || strcmp(disttype,'Uniform') ==1 || strcmp(disttype,'Bigaussian') ==1 
    figure
    yyaxis left
    stem(linspace(0,1,length(optweightfromAIC100)),optweightfromAIC100,'--square','LineWidth',2,'MarkerSize',10,'Color','#0072BD')
    hold on
    stem(linspace(0,1,length(optweightfromAIC25)),optweightfromAIC25*(r100to25),'--*','LineWidth',2,'MarkerSize',10,'Color','#6B9C28')
    hold on
    stem(linspace(0,1,length(optweightfromAIC10)),optweightfromAIC10*(r100to10),'--v','LineWidth',2,'MarkerSize',10,'Color','#E6AB1A')
    hold on
    stem(linspace(0,1,length(optweightfromAIC5)),optweightfromAIC5*(r100to5),'--x','LineWidth',2,'MarkerSize',10,'Color','k')
    hold on
    ylabel('Recovered Proportion of Population')
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[0 limsy(2)]);
    yyaxis right
    plot(linspace(0,1,101),sprobs100,'-','LineWidth',2,'MarkerSize',10,'Color','#A2142F')
    ylabel('Original Proportion of Population')
    legend('Recovered, N_t=100','Recovered, N_t=25','Recovered, N_t=10','Recovered, N_t=5','Original','Location','best','FontSize',12)
    set(gca,"FontSize",20)
    xlabel('Sensitivity to Treatment {\it s}')
end

PDFfiglabel=strcat(disttype,'N',noisestr,'_timepoints_PDF','.jpg');
saveas(gcf,PDFfiglabel);
PDFfiglabel=strcat(disttype,'N',noisestr,'_timepoints_PDF','.fig');
saveas(gcf,PDFfiglabel);


%plot scaled pmfs
if strcmp(disttype,'Normal') == 1 || strcmp(disttype,'Uniform') ==1 || strcmp(disttype,'Bigaussian') ==1 
    figure
    stem(linspace(0,1,length(optweightfromAIC100)),optweightfromAIC100,'--square','LineWidth',2,'MarkerSize',10,'Color','#0072BD')
    hold on
    stem(linspace(0,1,length(optweightfromAIC25)),optweightfromAIC25*(r100to25),'--*','LineWidth',2,'MarkerSize',10,'Color','#6B9C28')
    hold on
    stem(linspace(0,1,length(optweightfromAIC10)),optweightfromAIC10*(r100to10),'--v','LineWidth',2,'MarkerSize',10,'Color','#E6AB1A')
    hold on
    stem(linspace(0,1,length(optweightfromAIC5)),optweightfromAIC5*(r100to5),'--x','LineWidth',2,'MarkerSize',10,'Color','k')
    hold on
    ylabel('Recovered Proportion of Population')
    plot(linspace(0,1,101),sprobs100*length(sprobs100)/(sum(sprobs100)*length(optweightfromAIC100)),'-','LineWidth',2,'MarkerSize',10,'Color','#A2142F')
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[0 limsy(2)]);
    ylabel('Original Proportion of Population')
    legend('Recovered, N_t=100','Recovered, N_t=25','Recovered, N_t=10','Recovered, N_t=5','Original','Location','best','FontSize',12)
    set(gca,"FontSize",20)
    xlabel('Sensitivity to Treatment {\it s}')
    PDFfiglabel=strcat(disttype,'N',noisestr,'_timepoints_PDF_samescale','.jpg');
    saveas(gcf,PDFfiglabel);
    PDFfiglabel=strcat(disttype,'N',noisestr,'_timepoints_PDF_samescale','.fig');
    saveas(gcf,PDFfiglabel);
end


end