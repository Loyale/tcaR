function      [score,details]=contrast_tca_kde_multupdate(x,W,param,js,details);

if (nargout>1), detail=1; else detail=0; end

[m N]=size(x);
x=W*x;
sigma=x*x'/N;
Mut2ds=details.Mut2ds;
H2ds=details.H2ds;
weights=details.weights;
Hs=details.Hs;
ngrid=param.ngrid;
h=param.h;

lambdaC=param.lambdaC;
lambdaT=param.lambdaT;
lambdaG=param.lambdaG;

ks=1:m; ks(js)=[];
njs=length(js);
other=[];
for i=1:njs, other=[other js(i)*ones(1,m-njs)]; end
is=[repmat(ks,1,njs); other];
for i=2:njs, for j=1:i-1,    is = [is [ js(i); js(j)] ];   end, end

%reduces the set if we use a fixed tree
toremove=[]; 
if param.fixedtree==1,
    if isempty(param.tree)
        is=[];
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

notdone=ones(m,1);
notdone(ks)=0;

for k=1:size(is,2)
    
    [xs,ys,fs]=kdefft2d(x([is(1,k) is(2,k)],:),ngrid,h);
    fs=fs';  
    mu=(fs+1e-14).*log(fs+1e-14);
    H=-trapz(xs,trapz(ys,mu,2),1);
    H2ds(is(1,k),is(2,k))=H;
    H2ds(is(2,k),is(1,k))=H2ds(is(1,k),is(2,k));
    if notdone(is(2,k))
        fm=trapz(xs,fs,1);
        xm=ys;
        Hs(is(2,k))=trapz(xm,-(fm+1e-14).*log(fm+1e-14));
        notdone(is(2,k))=0;
    end
    if notdone(is(1,k))
        fm=trapz(xs,fs,2);
        xm=ys;
        Hs(is(1,k))=trapz(xm,-(fm+1e-14).*log(fm+1e-14));
        notdone(is(1,k))=0;
    end
    Mut2ds(is(1,k),is(2,k))=-H+Hs(is(1,k))+Hs(is(2,k));
    weights(is(1,k),is(2,k))=Mut2ds(is(1,k),is(2,k));
    corr=sigma(is(1,k),is(2,k))/sigma(is(1,k),is(1,k))^.5/sigma(is(2,k),is(2,k))^.5;
    weights(is(1,k),is(2,k))=weights(is(1,k),is(2,k))-lambdaG*.5*log(1+1e-14-corr^2)+ ...
        +lambdaC*.5*log(1+1e-14-corr^2);
    weights(is(2,k),is(1,k))=weights(is(1,k),is(2,k));
    Mut2ds(is(2,k),is(1,k))=Mut2ds(is(1,k),is(2,k));
end

for k=1:m
    if notdone(k),
        [xm,fm]=kdefft(x(k,:),ngrid,h);
        Hs(k)=trapz(xm,-(fm+1e-14).*log(fm+1e-14));
        notdone(k)=0;
    end
end

if param.fixedtree
        tree=param.tree;
    weightsum=0;
    for k=1:size(param.tree,2);
        weightsum = weightsum - weights(param.tree(1,k),param.tree(2,k));
    end
    score=sum(Hs)-log(abs(det(W)));
    score=score+lambdaT*weightsum;
else
if isfield(param,'edgeprior'), edgeprior=param.edgeprior; else edgeprior=0; end
[tree,weightsum] = minimum_nonspanning_tree_edges(edgeprior-weights,[],param.maxedges);
    score=sum(Hs)-log(abs(det(W)));
    score=score+lambdaT*weightsum;
end

if (detail),  details.Mut2ds=Mut2ds; details.weights=weights;  details.H2ds=H2ds; end
if (detail),  details.Hs=Hs; end
if (detail), details.tree=tree; end
