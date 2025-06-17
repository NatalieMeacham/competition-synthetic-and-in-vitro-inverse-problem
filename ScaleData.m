function [meanmatnorm,concvec,errmat]=ScaleData(choosedata)

%Scale data between 0 and 1

%INPUTS
    %choosedata: choice of dataset

%OUTPUTS
    %meanmatnorm: matrix of normalized data averaged across dosage
    %replicates
    %concvec: vector of concentration levels
    %errmat: matrix of error values for error bars 

load('MONOCLONAL_DATA.mat');
load('Mixture_Data.mat');

if choosedata=='R250'
    T=RESISTANT_250_BF;
elseif choosedata=='R500'
    T=RESISTANT_500_BF;
elseif choosedata=='S500'
    T=SENSITIVE_500_BF;
elseif choosedata=='S100'
    T=SENSITIVE_1000_BF;
elseif choosedata=='BF11'
    T=BF_11;
elseif choosedata=='BF12'
    T=BF_12;
elseif choosedata=='BF21'
    T=BF_21;
elseif choosedata=='BF41'
    T=BF_41;
else
    disp('??')
end

sz = size(T);
tsize=sz(3);
concsize=sz(2);
tvec=linspace(1,tsize,tsize);
concvec=linspace(1,concsize,concsize);
meanvec = zeros(tsize,1);
meanmat = zeros(concsize,tsize);
concnum=length(concvec);
sdmat = zeros(concsize,tsize);
errmat = zeros(concsize,tsize);

%find max value across monoclonal datasets and multiply it by 10 
maxval1=max(RESISTANT_250_BF);
maxval2=max(RESISTANT_500_BF);
maxval3=max(SENSITIVE_500_BF);
maxval4=max(SENSITIVE_1000_BF);
maxlist=[maxval1,maxval2, maxval3,maxval4];
maxval=10*max(maxlist,[],'all');

meanmatnorm=zeros(size(meanmat));

%normalize dataset using maxval
Tnormalized = T./maxval;

%compute normalized mean data
for a=concvec
    for t=tvec
        [s,m] = std(Tnormalized(:,a,t),'omitnan');
        meanmatnorm(a,t)= m;
        sdmat(a,t) = s;
        errmat(a,t) = s/sqrt(length(T(:,a,t)) - sum(isnan(T(:,a,t))));
    end 
end