%function[gls_optpar,converge_flag,AIC_GLS,weightedsol,rsgrid,finalerr]=GLSInverseFnC(rho, k, y0,tspan,rpoints,data)
%want data to be R500 dosage 1 mean

data=[0.0071 0.0081 0.0097 0.0121 0.0148 0.0184 0.0213 0.0264 0.0318 0.0388 0.0466 0.0534 0.0621 0.0723];
rho=0.0682
k=0.0077
tspan=linspace(9,48,14)
rpoints=4;
y0=data(1)

[gls_optpar,converge_flag,AIC_GLS,weightedsol,rsgrid,finalerr]=GLSInverseFnCData(rho, k, y0,tspan,rpoints,data)

% choosedata='S100'
% figures="yes"
% 
% [AICminvec, trymatoptrsgridS, optmatS, cdfRmatS, wsolmat,evvec,concvecS]=CompIPFn(choosedata,figures)