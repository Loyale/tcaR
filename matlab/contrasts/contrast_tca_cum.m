function [score,details]=contrat_tca_cum(x,W,param);

% ASSUMES THAT THE DATA HAS ZERO MEAN
% DOES NOT ASSUME UNIT VARIANCE

[m N]=size(x);
if (nargout>1), detail=1; else detail=0; end
x=W*x;
sigma=x*x'/N;


lambdaC=param.lambdaC;
lambdaG=param.lambdaG;
lambdaT=param.lambdaT;

notdone=ones(1,m);
Mut2ds=zeros(m);
weights=zeros(m);
H2ds=zeros(m);
H=zeros(1,m);

cum4=zeros(1,m);
cum4=sum(x.^4,2)/N;
cum3=sum(x.^3,2)/N;
cum2=sum(x.^2,2)/N;
Hs = -1 / 48 * ( cum4 ./ ( cum2 .^ 2 ) - 3 ) .^ 2 - 1 / 12 * cum3 .^ 2 ./ cum2 .^ 3 ; %+ 5/8 * cum3.^2 .* cum4 + 1/16 * (cum4-3).^3;

for i=2:m
    for j=1:i-1
        y=x([i j],:);
        sigy=sigma([i j],[i j]);
        y=inv(sqrtm(sigy))*y;
        Mut2ds(i,j)=-.5*log(1-sigma(i,j)^2 / sigma(i,i) / sigma(j,j) ) ...
            + 1/12*( sum( (y(1,:).^3).*y(2,:),2)/N )^2 ...
            + 1/8*(  sum( (y(1,:).^2).*(y(2,:).^2),2)/N -1)^2 ...
            + 1/12*( sum( (y(2,:).^3).*y(1,:),2)/N)^2 ...
            + 1/48*( (sum(y(1,:).^4)/N-3)^2 + (sum(y(2,:).^4)/N-3)^2) ...
            + 1/4*( sum( (y(2,:).^2).*y(1,:),2)/N)^2 ...
            + 1/4*( sum( (y(1,:).^2).*y(2,:),2)/N)^2 ...
            + 1/12*( (sum(y(1,:).^3)/N)^2 + (sum(y(2,:).^3)/N)^2) ...
            + Hs(i)+Hs(j);
        
        Mut2ds(j,i)=Mut2ds(i,j);
        weights(i,j)=Mut2ds(i,j);
        corr=sigma(i,j)/sigma(i,i)^.5/sigma(j,j)^.5;
        weights(i,j)=weights(i,j)+lambdaC*.5*log(1+1e-14-corr^2);
        weights(j,i)=weights(i,j);
    end
end

if isfield(param,'edgeprior'), edgeprior=param.edgeprior; else edgeprior=0; end
[tree,weightsum] = minimum_nonspanning_tree_edges(edgeprior-weights,[],param.maxedges);
score=sum(Hs)-log(abs(det(W)));
score=score+lambdaT*weightsum;
if (detail),  details.Mut2ds=Mut2ds; details.weights=weights;  end
if (detail),  details.Hs=Hs; end
if (detail), details.tree=tree; end
