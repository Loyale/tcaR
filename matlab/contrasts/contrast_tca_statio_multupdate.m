function      [score,details]=contrast_tca_statio_multupdate(x,W,param,js,details);

if (nargout>1), detail=1; else detail=0; end

[m N]=size(x);
x=W*x;
sigma=x*x'/N;
Mut2ds=details.Mut2ds;
weights=details.weights;
Hs=details.Hs;

nomega=size(param.periodo,3);
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

for k=1:length(js)
    i=js(k);
    for ko=1:nomega
        detspectral(ko)=real(W(i,:)*squeeze(param.periodo(:,:,ko))*W(i,:)');
        assert( detspectral(ko)>0,'negative spectral density');
    end
    Hs(i)=.5*sum(log(real(detspectral)))/(nomega);
end



for k=1:size(is,2)
        i=is(1,k); j=is(2,k);
    for ko=1:nomega
        detspectral(ko)=real(det(W([i j],:)*squeeze(param.periodo(:,:,ko))*W([i j],:)'));
        assert( detspectral(ko)>0,'negative spectral density');
    end
    Mut2ds(i,j)=Hs(i)+Hs(j)-.5*sum(log(real(detspectral)))/(nomega);
    Mut2ds(j,i)=Mut2ds(i,j);
    weights(i,j)=Mut2ds(i,j);
    corr=sigma(i,j)/sigma(i,i)^.5/sigma(j,j)^.5;     
    weights(i,j)=weights(i,j)+lambdaC*.5*log(1+1e-14-corr^2);
    weights(j,i)=weights(i,j);
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

if (detail),  details.Mut2ds=Mut2ds; details.weights=weights;  end
if (detail),  details.Hs=Hs; end
if (detail), details.tree=tree; end
