function [score,details]=contrat_tca_kde(x,W,param);
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

ngrid=param.ngrid;
h=param.h;
% have to compute all paiwise mutual informations
for i=2:m
    for j=1:i-1
        [xs,ys,fs]=kdefft2d(x([j i],:),ngrid,h);
        fs=fs';      % error of the kdefft function. TO CHECK!!
        mu=(fs+1e-14).*log(fs+1e-14);
        H=-trapz(xs,trapz(ys,mu,2),1);
        %[H,xs,ys,fs]=entropy2d_kdefft(x([j i],:),h,ngrid);
        H2ds(i,j)=H;
        H2ds(j,i)=H2ds(i,j);
        
        if (notdone(i))
            fm=trapz(xs,fs,1);
            xm=ys;
            Hs(i)=trapz(xm,-(fm+1e-14).*log(fm+1e-14));
            if (detail), details.margx{i}=xm; details.margf{i}=fm; end
            notdone(i)=0;
        end
        
        if (notdone(j))
            fm=trapz(ys,fs,2);
            xm=xs;
            Hs(j)=trapz(xm,-(fm+1e-14).*log(fm+1e-14));
            if (detail), details.margx{j}=xm; details.margf{j}=fm; end
            notdone(j)=0;
        end
        
        Mut2ds(i,j)=-H+Hs(i)+Hs(j);
        Mut2ds(j,i)=Mut2ds(i,j);
        weights(i,j)=Mut2ds(i,j);
        corr=sigma(i,j)/sigma(i,i)^.5/sigma(j,j)^.5;
        
        weights(i,j)=weights(i,j)+ lambdaC*.5*log(1+1e-14-corr^2);
        weights(j,i)=weights(i,j);
        
        if (detail)
            details.marg2dx{i,j}=xs;
            details.marg2dy{i,j}=ys;
            details.marg2df{i,j}=fs;
        end
    end
end
if isfield(param,'edgeprior'), edgeprior=param.edgeprior; else edgeprior=0; end
[tree,weightsum] = minimum_nonspanning_tree_edges(edgeprior-weights,[],param.maxedges);
score=sum(Hs)-log(abs(det(W)));
score=score+lambdaT*weightsum;

if (detail),  details.Mut2ds=Mut2ds; details.weights=weights;  details.H2ds=H2ds; end
if (detail),  details.Hs=Hs; end
if (detail), details.tree=tree; end
