sgrid=linspace(0,1,101);
sprobs=DistFn2('OnePoint',sgrid,0,1);
disttype='Other';
ttreatment=linspace(0,10,11);
tgrowth=linspace(0,5,6);
y0=0.1;
treatcycles=9;
growthcycles=8; %growth cycles must be either same as treatcycles or one fewer
noisesize=0.0;
rho=0.3;
k=0.45;
[tmat,gmat,ttreatmat,tgrowthmat,sprobsmat] = SyntheticCyclesFn(sgrid,sprobs,ttreatment,tgrowth,y0,treatcycles,growthcycles,noisesize,rho,k);


% figure
% for i=1:treatcycles
%     plot(sgrid,sprobsmat(:,i)./sum(sprobsmat(:,i)))
%     hold on
% end
% rpoints=7;
% [gls_optpar,converge_flag,AIC_GLS,weightedsol,rsgrid]=GLSInverseFnCCycles(rho,k,y0,ttreatment,rpoints,tmat(:,1)')
% % %run IP on each set of treatment data
% [sgrid,sprobs, rsgrid,optweightfromAIC,t,propdata,sweightedsol,cdfS,cdfR]= GLSInverseScriptC(tmat(:,1),ttreatmat(:,1),'y')
% 
% figure
% rsgrid=linspace(0,1,rpoints);
% plot(sgrid,sprobs)
% hold on
% plot(rsgrid,gls_optpar)

propdata=tmat(:,1)
tspan=ttreatment;
figures='y';
[sgrid,sprobs, rsgrid,optweightfromAIC,t,propdata,sweightedsol,cdfS,cdfR]= GLSInverseScriptCCycles(sgrid,sprobs,rho,k,y0,propdata',tspan,'y',disttype)


%run IP on each treatment cycle
optweightmat=zeros(30,treatcycles); %matrix of recovered values
rsgridvec=zeros(1,treatcycles); %number of nodes in each recovered dist
sweightedsolmat=zeros(length(tspan),treatcycles); %mat w sweightedsol for each cycle
cdfRmat=zeros(30,treatcycles); %matrix of recovered cdfS
for i=1:treatcycles
    [~,~, rsgrid,optweightfromAIC,t,propdata,sweightedsol,cdfS,cdfR]= GLSInverseScriptCCycles(sgrid,sprobsmat(:,i),rho,k,tmat(1,i),tmat(:,i)',tspan,'y',disttype)
    optweightmat(1:length(rsgrid),i)=optweightfromAIC;
    rsgridvec(i)=length(rsgrid);
    sweightedsolmat(:,i)=sweightedsol;
    cdfRmat(1:length(rsgrid),i)=cdfR;
end
%% 

clr=autumn(treatcycles)
%'Color',clr(a,:)

figure
for i=1:treatcycles
    subplot(treatcycles,1,i)
    plot(sgrid,sprobsmat(:,i)./sum(sprobsmat(:,i)),'*-','Color',clr(i,:),'LineWidth',2)
    hold on
    plot(linspace(0,1,rsgridvec(i)),optweightmat(1:rsgridvec(i)),'--','Color',clr(i,:),'LineWidth',2)
    legendstring=strcat('c',num2str(i),'original')
    legend(legendstring,'recovered')
    %title(['Treatment Cycle', num2str(i) ])
end
%legend('original 1','recovered 1','original 2','recovered 2','original 3','recovered 3','original 4','recovered 4')

figure
legendstring2=[]
for i=1:treatcycles
    plot(sgrid,sprobsmat(:,i)./sum(sprobsmat(:,i)),'*-','Color',clr(i,:),'LineWidth',2)
    hold on
    legendstring2a=strcat('Cycle',{' '},num2str(i))
    legendstring2=[legendstring2,legendstring2a]
    %legend(legendstring,'recovered')
    %title(['Treatment Cycle', num2str(i) ])
end
legend(legendstring2)

figure
for i=1:treatcycles
    plot( ttreatmat(:,i),tmat(:,i),"k*")
    hold on
    plot(ttreatmat(:,i),sweightedsolmat(:,i),'k-')
    if i<=growthcycles
        plot(tgrowthmat(:,i), gmat(:,i),'r.','MarkerSize',12 )
        hold on
    end
end

% figure
% for i=1:treatcycles
%     % plot(sgrid,sprobsmat(:,i)./sum(sprobsmat(:,i)),'--','Color',clr(i,:))
%     % hold on
%     plot(linspace(0,1,rsgridvec(i)),optweightmat(1:rsgridvec(i)),'Color',clr(i,:))
%     hold on
% end

% figure
% plot(sgrid,sprobsmat(:,1)./sum(sprobsmat(:,1)),'--','Color',clr(1,:))
% hold on
% plot(linspace(0,1,rsgridvec(1)),optweightmat(1:rsgridvec(1)),'Color',clr(1,:))
% 
% hold on
% plot(sgrid,sprobsmat(:,2)./sum(sprobsmat(:,2)),'--','Color',clr(2,:))
% hold on
% plot(linspace(0,1,rsgridvec(2)),optweightmat(1:rsgridvec(2)),'Color',clr(2,:))
% 
% hold on
% plot(sgrid,sprobsmat(:,3)./sum(sprobsmat(:,3)),'--','Color',clr(3,:))
% hold on
% plot(linspace(0,1,rsgridvec(3)),optweightmat(1:rsgridvec(3)),'Color',clr(3,:))
% 
% hold on
% plot(sgrid,sprobsmat(:,4)./sum(sprobsmat(:,4)),'--','Color',clr(4,:))
% hold on
% plot(linspace(0,1,rsgridvec(4)),optweightmat(1:rsgridvec(4)),'Color',clr(4,:))
% figure
% plot(sgrid,sprobs)
% hold on
% plot(rsgrid,optweightfromAIC)
% 
% figure
% plot(tspan,propdata,"*")
% hold on
% plot(t,sweightedsol)