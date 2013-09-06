function      [score,details]=contrast_tca_cum_multupdate(x,W,param,js,details)
if (nargout>1), detail=1; else detail=0; end
[m N]=size(x);
x=W*x;
sigma=x*x'/N;
Mut2ds=details.Mut2ds;
weights=details.weights;
Hs=details.Hs;
lambdaC=param.lambdaC;
lambdaG=param.lambdaG;
lambdaT=param.lambdaT;


cum4=zeros(1,m);
cum4=sum(x.^4,2)/N;
cum3=sum(x.^3,2)/N;
Hs=-1/48 * (cum4 - 3).^2 - 1/12 * cum3.^2; % + 5/8 * cum3.^2 .* cum4 + 1/16 * (cum4 - 3).^3;

% is=1:m; is([i j])=[];
% is=[is is i ; i*ones(1,m-2) j*ones(1,m-1)];

ks=1:m; ks(js)=[];
njs=length(js);
other=[];
for i=1:njs, other=[other js(i)*ones(1,m-njs)]; end
is=[repmat(ks,1,njs); other];
for i=2:njs, for j=1:i-1,    is = [is [ js(i); js(j)] ];   end, end
%reduces the set if we use a fixed tree
toremove=[]; 
if param.fixedtree==1,
    if isempty(param.tree), is=[]; 
    else
    alltree=[ param.tree(1,:)+m*param.tree(2,:) param.tree(2,:)+m*param.tree(1,:) ];
    for k=1:size(is,2);
        if ~ismember(is(1,k)+m*is(2,k),alltree)
            toremove=[toremove k];
        end
    end
    is(:,toremove)=[];
end
end


for k=1:size(is,2)
    y=x([is(1,k) is(2,k)],:);
    sigy=sigma([is(1,k) is(2,k)],[is(1,k) is(2,k)]);
    y=inv(sqrtm(sigy))*y;
    Mut2ds(is(1,k),is(2,k))=-.5*log(1-sigma(is(1,k),is(2,k))^2 / sigma(is(1,k),is(1,k)) / sigma(is(2,k),is(2,k)) ) ...
        + 1/12*( sum( (y(1,:).^3).*y(2,:),2)/N )^2 ...
        + 1/8*(  sum( (y(1,:).^2).*(y(2,:).^2),2)/N -1)^2 ...
        + 1/12*( sum( (y(2,:).^3).*y(1,:),2)/N)^2 ...
        + 1/48*( (sum(y(1,:).^4)/N-3)^2 + (sum(y(2,:).^4)/N-3)^2) ...
        + 1/4*( sum( (y(2,:).^2).*y(1,:),2)/N)^2 ...
        + 1/4*( sum( (y(1,:).^2).*y(2,:),2)/N)^2 ...
        + 1/12*( (sum(y(1,:).^3)/N)^2 + (sum(y(2,:).^3)/N)^2) ...
        + Hs(is(1,k))+Hs(is(2,k));         
    weights(is(1,k),is(2,k))=Mut2ds(is(1,k),is(2,k));
    corr=sigma(is(1,k),is(2,k))/sigma(is(1,k),is(1,k))^.5/sigma(is(2,k),is(2,k))^.5;
    weights(is(1,k),is(2,k))=weights(is(1,k),is(2,k))-lambdaG*.5*log(1+1e-14-corr^2)+ ...
        +lambdaC*.5*log(1+1e-14-corr^2);
    weights(is(2,k),is(1,k))=weights(is(1,k),is(2,k));
    Mut2ds(is(2,k),is(1,k))=Mut2ds(is(1,k),is(2,k));
end


if param.fixedtree
    tree=param.tree;
    weightsum=0;
    for k=1:size(param.tree,2);
        weightsum = weightsum - weights(param.tree(1,k),param.tree(2,k));
    end
score=sum(Hs)-(1+lambdaG)*log(abs(det(W)));
    score=score+lambdaT*weightsum;
else
if isfield(param,'edgeprior'), edgeprior=param.edgeprior; else edgeprior=0; end
[tree,weightsum] = minimum_nonspanning_tree_edges(edgeprior-weights,[],param.maxedges);
score=sum(Hs)-(1+lambdaG)*log(abs(det(W)));
    score=score+lambdaT*weightsum;
end

if (detail),  details.Mut2ds=Mut2ds; details.weights=weights;  end
if (detail),  details.Hs=Hs; end
if (detail), details.tree=tree; end

