function [t, cmat,weightedsol] = ForwardFnFindRhoC(rho, y0, tspan_data)

tinit=tspan_data(1);
tfinal=tspan_data(end);
tspan_data;
tspan=linspace(tinit,tfinal,131);
h=tspan(2) - tspan(1);
sprobs=1;

%list of indices of tspan that also have a t_data point
tdindices = [];
for i=1:length(tspan)
    if ismember(tspan(i),tspan_data) == 1
        tdindices=[tdindices i];
    end
end

%tdindices
dcdt = @(csol,u) rho.*csol.*(1 - u);

tpoints=length(tspan);
points=length(sprobs);

uvec=zeros(tpoints,1);

y0vec=(y0*ones(points,1))';

uvec(1)=sum((sprobs/sum(sprobs)).*y0vec);

%uvec

%cmat=zeros(length(sgrid),length(tspan));
cmat=zeros(1,length(tspan));

disp(size(cmat));

    cmat(1,1)=y0;
    
    %cmat

for j=2:tpoints

    u=uvec(j-1);

    k1=h.*dcdt(cmat(:,j-1),u);
    k2=h.*dcdt(cmat(:,j-1),u + k1/2);
    k3=h.*dcdt(cmat(:,j-1),u + k2/2);
    k4=h.*dcdt(cmat(:,j-1),u + k3);
    cmat(:,j)=cmat(:,j-1) + (k1 + 2*k2 + 2*k3 + k4)/6;
   
    uvec(j)=((sprobs/sum(sprobs))*cmat(:,j));

end

cmat=cmat';
weightedsol=(sprobs./(sum(sprobs)))*cmat'; %old

t=tspan;
t=t';
uvec';

t=t([tdindices]);
cmat=cmat([tdindices],:);
weightedsol=weightedsol([tdindices]);
end