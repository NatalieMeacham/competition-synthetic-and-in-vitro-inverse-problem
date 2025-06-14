function [sgrid,sprobs, rsgrid,optweightfromAIC,t,propdata,sweightedsol,cdfS,cdfR]= GLSInverseScriptC(disttype, tpoints,noisesize,figures)
    
    %Declare variables 
    points=101; 
    tfinal=50;
    tspan=linspace(0,tfinal,tpoints);
    pointsstr=string(points);
    noisestr=string(noisesize);

    %Parameters for synthetic data
    rho=0.3;
    k=0.45;
    y0=.2;
    
    %Define sensitivity mesh
    a=0;
    b=1;
    sgrid=linspace(a,b,points);
    
    %Create original distribution, synthetic data, add noise
    [sprobs] = DistFn2(disttype,sgrid,a,b);
 
    %Generate synthetic data with proportional error
    [t, cmat,weightedsol] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan);
    propdata=weightedsol.*(1 + noisesize*randn(size(weightedsol)));

    cmatbeforeIP=cmat
    
    %%
    
    %Loop through IP over different meshes -> find optimal mesh
    rpointsvec=4:1:30;
    AICvec=zeros(length(rpointsvec),1);
    for i=1:length(rpointsvec)
        [gls_optpar,~,AIC_GLS,~,~,~,rsgrid]=GLSInverseFnC(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpointsvec(i),propdata);
        AICvec(i)=AIC_GLS;
    end
    
    [M,I]=min(AICvec);
    optrpoints=rpointsvec(I);
    
    %get optweight from OLSInverseFn for optrpoints
    [optweightfromAIC,~,~,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFnC(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,optrpoints,propdata);
    optweightfromAIC
    rsgrid
    %keyboard
    [t,cmat,sweightedsol] = RK4FunctionC(rsgrid, optweightfromAIC', rho, k, y0, tspan);
    sweightedsol;
    cmatafterIP=cmat
    size(cmatafterIP)
    cmatbeforeIP 
    size(cmatbeforeIP)

    % figure
    % for i=1:length(rsgrid)
    %     plot(t,cmatafterIP(:,i),'b')
    %     hold on
    % end
    % for j=1:length(sgrid)
    %     plot(t,cmatbeforeIP(:,j),'r')
    %     hold on
    % end

     %make and plot cdfs to accommodate cts dists:
    cdfR=zeros(optrpoints,1);
    cdfS=zeros(points,1);
    if strcmp(disttype,'Normal') == 1 || strcmp(disttype,'Uniform') ==1 || strcmp(disttype,'Bigaussian') ==1
        partialsumS=0;
        cdfS(1)=(1/points)*sprobs(1);
        for i=2:points
            partialsumS=trapz(sgrid(1:i),sprobs(1:i));
            cdfS(i)=partialsumS;
        end
        partialsumR=0;
        for i=1:optrpoints
            partialsumR=sum(optweightfromAIC(1:i));
            cdfR(i)=partialsumR;
        end
        disp('Here is a cdf for a continuous dist.')
    elseif strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1
        partialsumS=0;
        for i=1:points
            partialsumS=sum(sprobs(1:i));
            cdfS(i)=partialsumS;
        end
        partialsumR=0;
        for i=1:optrpoints
            partialsumR=sum(optweightfromAIC(1:i));
            cdfR(i)=partialsumR;
        end
        disp('Here is a cdf for a discrete dist.')
    end
    

if figures == 'y'

    figure
    plot(rpointsvec,AICvec,'LineWidth',2,'Color','blue')
    xlabel('Number of Points in Recovered Dist.')
    ylabel('AIC Score')
    %title('AIC Score for Differing Meshes')
    set(gca,"FontSize",20)
    hold on
    plot(rpointsvec(I),M,'m*','LineWidth',2,'MarkerSize',20)
    hold off
    legend('AIC','Min. AIC','Location','northeast')
    %AICfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'AIC','.jpg');
    %saveas(gcf,AICfiglabel);
    %AICfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'AIC','.fig');
    %saveas(gcf,AICfiglabel);
    
    
    %compare output of fmincon with original dist (no lines between recovered pts)
    figure
    yyaxis left
    stem(rsgrid,optweightfromAIC,'--b','MarkerSize',8,'LineWidth',2) %,'MarkerSize,'8,'LineWidth',2,'Color,''blue')
    ylabel('Recovered Proportion of Population')
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[0 limsy(2)]);
    ylimval=max(max(sprobs),max(optweightfromAIC));
    if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1
        ylim([0 ylimval])
    end
    if strcmp(disttype,'Uniform') == 1
        ylim([0 2.5*median(optweightfromAIC)])
    end
    hold on
    yyaxis right
    plot(sgrid,sprobs,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
    legend('Recovered','Original','Location','northeast')
    xlabel('Sensitivity to Treatment {\it s}')
    ylabel('Original Proportion of Population')
    if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1
        ylim([0 ylimval])
    end
    if strcmp(disttype,'Uniform') == 1
        ylim([0 2.5*median(sprobs)])
    end
    set(gca,"FontSize",20)
    %PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Dists','.jpg');
    %saveas(gcf,PMFfiglabel);
    %PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Dists','.fig');
    %saveas(gcf,PMFfiglabel);
    
   
    
    figure
    X=rsgrid;
    Y=cdfR;
    stairs(X,Y,'--o','MarkerSize',8,'LineWidth',2,'Color','blue')
    ylim([0 1])
    hold on
    plot(sgrid, cdfS,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
    xlabel('Sensitivity to Treatment {\it s}')
    ylabel('Cumulative Recovered Proportion')
    legend('Recovered','Original','Location','Southeast')
    set(gca,"FontSize",20)
    %CDFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'CDF','.jpg');
    %saveas(gcf,CDFfiglabel);
    %CDFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'CDF','.fig');
    %saveas(gcf,CDFfiglabel);
    
    
    %[t,~,sweightedsol] = RK4FunctionC(rsgrid, optweightfromAIC', rho, k, y0, tspan);
    %sweightedsol
    %sweightedsolGuess
    figure
    plot(t,sweightedsol,'LineWidth',2,'Color','blue')
    hold on
    plot(t,propdata','rd','LineWidth', 2,'MarkerSize',8)
    hold off
    legend('Est. Tumor Volume','Simulated Data','Location','southeast')
    xlabel('Time')
    ylabel('Aggregated Tumor Volume')
    %title('Weighted Sols for Original vs. Recovered Dists.')
    set(gca,"FontSize",20)
    %Fitfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Fit','.jpg');
    %saveas(gcf,Fitfiglabel);
    %Fitfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Fit','.fig');
    %saveas(gcf,Fitfiglabel);

    %keyboard
    
    
    %compare output of fmincon with original dist on same scaled axes (no lines between recovered pts)
    if strcmp(disttype,'Uniform') == 1 || strcmp(disttype,'Normal') == 1 || strcmp(disttype,'Bigaussian') == 1
        figure
        stem(rsgrid,optweightfromAIC,'--b','MarkerSize',8,'LineWidth',2) %,'MarkerSize,'8,'LineWidth',2,'Color,''blue')
        ylabel('Proportion of Population')
        hold on
        plot(sgrid,sprobs*length(sprobs)/(sum(sprobs)*length(optweightfromAIC)),'--o','MarkerSize',8,'LineWidth',2,'Color','red')
        limsy=get(gca,'YLim');
        set(gca,'Ylim',[0 limsy(2)]);
        legend('Recovered','Original','Location','northeast')
        xlabel('Sensitivity to Treatment {\it s}')
        set(gca,"FontSize",20)
        %PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Dists_samescale','.jpg');
        %saveas(gcf,PMFfiglabel);
        %PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Dists_samescale','.fig');
        %saveas(gcf,PMFfiglabel);
    end

    if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') == 1 
        figure
        stem(rsgrid,optweightfromAIC,'--b','MarkerSize',8,'LineWidth',2) 
        limsy=get(gca,'YLim');
        set(gca,'Ylim',[0 limsy(2)]);
        hold on
        plot(sgrid,sprobs,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
        legend('Recovered','Original','Location','northeast')
        xlabel('Sensitivity to Treatment {\it s}')
        ylabel('Proportion of Population')
        ylimval=max(max(sprobs),max(optweightfromAIC));
        ylim([0 ylimval])
        set(gca,"FontSize",20)
        %PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'DistsOneAxis','.jpg');
        %saveas(gcf,PMFfiglabel);
        %PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'DistsOneAxis','.fig');
        %saveas(gcf,PMFfiglabel);
    end

end
end
