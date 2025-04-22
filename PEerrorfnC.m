function[err,bestapprox]=PEerrorfn(s, gencmatC, data,weights)

bestapprox=(gencmatC*s)'; %bestapprox gets weightedsol by multiplying s*csol

%err=sum(weights.*((bestapprox - constdata).^2));
err=sum((weights.^(1)).*((bestapprox - data)).^2);
end