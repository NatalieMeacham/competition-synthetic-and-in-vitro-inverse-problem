function [meanmatnorm,isnanmat, tspanvec, tspanmat,concvecA]=DropNaNsFn(concvec,meanmatnorm,tspan)

%Remove NaN values from the data and if a concentration of data has too
%many NaNs, remove that row and adjust the resulting mean data matrix and
%vector of concentration options 

%INPUTS
    %concvec: vector of concentration numbers 
    %meanmatnorm: matrix of normalized mean data 
    %tspan: time vector

%OUTPUTS
    %meanmatnorm: adjusted maxtrix of normalized mean data
    %isnanmat: matrix describing NaN values
    %tspanvec: vector describing time spans for each concentration
    %tspanmat: matrix describing time spans for each concentration
    %concvecA: adjusted set of concentration numbers 

isnanmat=zeros(size(meanmatnorm));
tspanvec=zeros(length(concvec),1);
tspanmat=zeros(length(concvec),length(tspan));

for a=concvec
    tspanmat(a,:)=tspan;
end

for a=concvec
    isnanmat(a,:)=isnan(meanmatnorm(a,:));
    k = find(isnanmat(a,:));
    if isempty(k)==1
        tspanvec(a)=14;
    else
        tfillvec=tspanmat(a,:);
        tfillvec(k)=[];
        tspanmat(a,1:14-length(k))=tfillvec;
        tspanmat(a,15-length(k):end)=zeros;
        tspanvec(a)=14-length(k);
    end
    B = rmmissing(meanmatnorm(a,:));
    meanmatnorm(a,1:length(B))=B;
    meanmatnorm(a,(length(B)+1):end)=zeros(1,length(tspan)-length(B));
end

%for data missing too many NaNs, remove that row and adjust all dimensions
for a=concvec
    if sum(isnanmat(a,:))>length(isnanmat(a,:))/2
        meanmatnorm(a,:)=[];
        tspanvec(a)=[];
        tspanmat(a,:)=[];
        concvec(a)=[];
        concvecA=concvec;
    else
        concvecA=concvec;
    end
end