function [sprobs] = DistFn2(disttype,sgrid,a,b)

%FUTURE UPDATE: do "BigaussianLow" and "BigaussianHigh" for two
%different bigaussian dists

if string(disttype) == 'Normal' %this should be cts
    mu = .5; 
    sigma = .09;
    sprobs = pdf('Normal',sgrid,mu,sigma); 
    sprobs = sprobs/(trapz(sgrid,sprobs));
% elseif string(disttype) == 'CtsNormal'
%     mu = .5; 
%     sigma = .09;
%     %sprobs = pdf('Normal',sgrid,mu,sigma); 
%     %normal = @(x) (1/(sigma*sqrt(2*pi)))*exp(-.5*((x - mu)/sigma).^2);
%     %sprobs=normal(sgrid)
%     sprobs = normpdf(sgrid,mu,sigma); 
 elseif string(disttype) == 'Bigaussian' %this should be cts 
     % %main paper figure (equally weighted bumps with different SDs)
     % mu1=.3;
     % mu2=.7; 
     % sigma1=.05;
     % sigma2=.1; %originally 0.1
     % weight1 = 0.5; %paper new version
     % weight2 = 0.5;

     % %original figure (wider, more heavily weighted left bump)
     mu1=.3;
     mu2=.7; 
     sigma1=.1;
     sigma2=.05;
     weight1=.9; %typically use this pair
     weight2=.1;

     % % % % %original figure with sides swapped
     % mu1=.3;
     % mu2=.7; 
     % sigma2=.1;
     % sigma1=.05;
     % weight2=.9; %typically use this pair
     % weight1=.1;

     % %equally sized and weighted bumps
     % mu1=.3;
     % mu2=.7; 
     % sigma1=.05;
     % sigma2=.05;
     % weight1=.5; %typically use this pair
     % weight2=.5;

     % %narrower equally sized and weighted bumps
     % mu1=.3;
     % mu2=.7; 
     % sigma1=.01;
     % sigma2=.01;
     % weight1=.5; %typically use this pair
     % weight2=.5;

     % %equally sized and weighted bumps farther apart 
     % mu1=.15;
     % mu2=.85; 
     % sigma1=.05;
     % sigma2=.05;
     % weight1=.5; %typically use this pair
     % weight2=.5;

     bigaussian = @(x) weight1*(1/(sigma1*sqrt(2*pi)))*exp(-.5*((x - mu1)/sigma1).^2) + weight2*(1/(sigma2*sqrt(2*pi)))*exp(-.5*((x - mu2)/sigma2).^2);
     sprobs=bigaussian(sgrid);
     sprobs = sprobs/(trapz(sgrid,sprobs));        % Truncated to ensure it is a odf
elseif string(disttype) == 'Uniform'
     sprobs = pdf('Uniform',sgrid, a, b);
     sprobs = sprobs/(trapz(sgrid,sprobs));
elseif string(disttype) == 'OnePoint'
    pointval=0.4; %usually 0.4
    %point=pointval*(length(sgrid) - 1);
    %point=(sgrid(pointval*100 + 1));
    sprobs=zeros(length(sgrid),1);
    sprobs(pointval*100)=1;%Note to self: the pointval*100 is the problem for different numbers of points
    sprobs=sprobs';
    sprobs = sprobs/(sum(sprobs));
elseif string(disttype) == 'TwoPoints'
    pointval1=0.3; %normally .3
    pointval2=0.8;
    prop1=0.9; %normally 0.9
    prop2=1 - prop1;
    sprobs=zeros(length(sgrid),1);
    sprobs(pointval1*100)=prop1;
    sprobs(pointval2*100)=prop2;
    sprobs=sprobs';
    sprobs = sprobs/(sum(sprobs));
else 
    disp('Choose a distribution that works.')
end

%Normalize sprobs to sum to 1 for whichever distribution:

end
