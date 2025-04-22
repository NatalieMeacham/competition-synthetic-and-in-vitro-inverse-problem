function [sprobs] = DistFn(disttype,sgrid,a,b)

%FUTURE UPDATE: do "BigaussianLow" and "BigaussianHigh" for two
%different bigaussian dists

if string(disttype) == 'Normal'
    mu = .5; 
    sigma = .09;
    sprobs = pdf('Normal',sgrid,mu,sigma); 
 elseif string(disttype) == 'Bigaussian'
     mu1=.3;
     mu2=.7; 
     % sigma1=.05;
     % sigma2=.1;
     % weight1=.1;
     % weight2=.9;
     sigma1=.1;
     sigma2=.05;
     weight1=.9;
     weight2=.1;
     % weight1=.1;
     % weight2=.9;
     %  %bigaussian = @(x) (1/(sigma1*sqrt(2*pi)))*exp(-.5*((x - mu1)/sigma1).^2) + (1/(sigma2*sqrt(2*pi)))*exp(-.5*((x - mu2)/sigma2).^2);
     bigaussian = @(x) weight1*(1/(sigma1*sqrt(2*pi)))*exp(-.5*((x - mu1)/sigma1).^2) + weight2*(1/(sigma2*sqrt(2*pi)))*exp(-.5*((x - mu2)/sigma2).^2);
     sprobs=bigaussian(sgrid);
elseif string(disttype) == 'Uniform'
     sprobs = pdf('Uniform',sgrid, a, b);
elseif string(disttype) == 'OnePoint'
    pointval=0.4;
    %point=pointval*(length(sgrid) - 1);
    %point=(sgrid(pointval*100 + 1));
    sprobs=zeros(length(sgrid),1);
    sprobs(pointval*100)=1;
    sprobs=sprobs';
elseif string(disttype) == 'TwoPoints'
    pointval1=0.3;
    pointval2=0.8;
    prop1=0.4;
    prop2=1 - prop1;
    sprobs=zeros(length(sgrid),1);
    sprobs(pointval1*100)=prop1;
    sprobs(pointval2*100)=prop2;
    sprobs=sprobs';
else 
    disp('Choose a distribution that works.')
end

%Normalize sprobs to sum to 1 for whichever distribution:
sprobs = sprobs/(sum(sprobs));
end