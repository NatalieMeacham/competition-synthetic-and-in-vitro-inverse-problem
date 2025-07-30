%% run inverse problem for one chosen dataset
%change dataset to get figures 5-7 from main paper 
choosedata=['BF21']
figures = ['y']
[AICminvec, trymatoptrsgridS, optmatS, cdfRmatS, wsolmat,evvec,concvecS]=GLSAllDataInclNaNsC(choosedata,figures);

%% run inverse problem for all datasets (no figures)
%%R250, R500, S500, and S1000 or BF11,BF12,BF21,BF41
figures = ['n'];
[AICminvecR250, trymatoptrsgridSR250, optmatSR250, cdfRmatSR250, wsolmatR250,evvecR250,concvecSR250]=GLSAllDataInclNaNsC('R250',figures)
[AICminvecR500, trymatoptrsgridSR500, optmatSR500, cdfRmatSR500, wsolmatR500,evvecR500,concvecSR500]=GLSAllDataInclNaNsC('R500',figures)
[AICminvecS500, trymatoptrsgridSS500, optmatSS500, cdfRmatSS500, wsolmatS500,evvecS500,concvecSS500]=GLSAllDataInclNaNsC('S500',figures)
[AICminvecS1000, trymatoptrsgridSS1000, optmatSS1000, cdfRmatSS1000, wsolmatS1000,evvecS1000,concvecSS1000]=GLSAllDataInclNaNsC('S100',figures)
[AICminvecBF11, trymatoptrsgridSBF11, optmatSBF11, cdfRmatSBF11, wsolmatBF11,evvecBF11,concvecSBF11]=GLSAllDataInclNaNsC('BF11',figures)
[AICminvecBF12, trymatoptrsgridSBF12, optmatSBF12, cdfRmatSBF12, wsolmatBF12,evvecBF12,concvecSBF12]=GLSAllDataInclNaNsC('BF12',figures)
[AICminvecBF21, trymatoptrsgridSBF21, optmatSBF21, cdfRmatSBF21, wsolmatBF21,evvecBF21,concvecSBF21]=GLSAllDataInclNaNsC('BF21',figures)
[AICminvecBF41, trymatoptrsgridSBF41, optmatSBF41, cdfRmatSBF41, wsolmatBF41,evvecBF41,concvecSBF41]=GLSAllDataInclNaNsC('BF41',figures)

%% plot mean sensitivity for all 8 dists (figure 8 from main paper)
dosevec= [0.01, 0.0260,0.0536, 0.1094, 0.2244, 0.3410,0.4603,1.1908,2.4430,3.7133,5.0119];

clr=jet(6);
figure
semilogx(dosevec(concvecSR250),evvecR250,'d','MarkerSize',8,'LineWidth',2,'Color',clr(1,:),'MarkerFaceColor',clr(1,:))
hold on
semilogx(dosevec(concvecSR500),evvecR500,'o','MarkerSize',8,'LineWidth',2,'Color',clr(1,:),'MarkerFaceColor',clr(1,:))
hold on
semilogx(dosevec(concvecSBF12),evvecBF12,'d','MarkerSize',8,'LineWidth',2,'Color',clr(2,:),'MarkerFaceColor',clr(2,:))
hold on
semilogx(dosevec(concvecSBF11),evvecBF11,'d','MarkerSize',8,'LineWidth',2,'Color',clr(3,:),'MarkerFaceColor',clr(3,:))
hold on
semilogx(dosevec(concvecSBF21),evvecBF21,'d','MarkerSize',8,'LineWidth',2,'Color',clr(4,:),'MarkerFaceColor',clr(4,:))
hold on
semilogx(dosevec(concvecSBF41),evvecBF41,'d','MarkerSize',8,'LineWidth',2,'Color',clr(5,:),'MarkerFaceColor',clr(5,:))
hold on
semilogx(dosevec(concvecSS500),evvecS500,'d','MarkerSize',8,'LineWidth',2,'Color',clr(6,:),'MarkerFaceColor',clr(6,:))
hold on
semilogx(dosevec(concvecSS1000),evvecS1000,'o','MarkerSize',8,'LineWidth',2,'Color',clr(6,:),'MarkerFaceColor',clr(6,:))
legend('Resistant 1','Resistant 2','S:2R Mixture','S:R Mixture','2S:R Mixture','4S:R Mixture','Sensitive 1','Sensitive 2','Location','northwest','FontSize',16)
xlabel('Log(Dosage)')
ylabel('Mean Recovered Sensitivity')
ylim([0 1.02])
set(gca,"FontSize",20)
saveas(gcf,'Figures/MeanRecoveredSensitivity.fig');
saveas(gcf,'Figures/MeanRecoveredSensitivity.png');