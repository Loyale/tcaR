function [score,details]=contrat_tca_statio(x,W,param);

% ASSUMES THAT THE DATA HAS ZERO MEAN

[m N]=size(x);
if (nargout>1), detail=1; else detail=0; end
x=W*x;
sigma=x*x'/N;


nomega=size(param.periodo,3);
weights=zeros(m);
Hs=zeros(1,m);
for i=1:m
    for k=1:nomega
        detspectral(k)=real(W(i,:)*squeeze(param.periodo(:,:,k))*W(i,:)');
        assert( detspectral(k)>0,'negative spectral density');
    end
    Hs(i)=.5*sum(log(real(detspectral)))/(nomega);
end


lambdaC=param.lambdaC;
lambdaG=param.lambdaG;
lambdaT=param.lambdaT;
for i=2:m
    for j=1:i-1
        for k=1:nomega
            detspectral(k)=real(det(W([i j],:)*squeeze(param.periodo(:,:,k))*W([i j],:)'));
            assert( detspectral(k)>0,'negative spectral density');
        end
        Mut2ds(i,j)=Hs(i)+Hs(j)-.5*sum(log(real(detspectral)))/(nomega);
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
