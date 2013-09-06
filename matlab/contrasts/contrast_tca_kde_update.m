function      [score,details]=contrast_tca_kde_update(x,W,param,i,details);

if (param.fixedtree)
   fixedtree=1;
   tree=param.tree;
else fixedtree=0; end

[m N]=size(x);
x=W'*x;
sigma=x*x'/N;
Mut2ds=details.Mut2ds;
H2ds=details.H2ds;
Hs=details.Hs;
ngrid=param.ngrid;
h=param.h;

lambdaC=param.lambdaC;
lambdaG=param.lambdaG;


if (~fixedtree)
   notdone=zeros(1,m);
   notdone=1;
   for k=1:m
      if (k~=i)
         [xs,ys,fs]=kdefft2d(x([k i],:),ngrid,h);
         fs=fs';  
         mu=(fs+1e-14).*log(fs+1e-14);
         H=-trapz(xs,trapz(ys,mu,2),1);
         H2ds(i,k)=H;
         H2ds(k,i)=H2ds(i,k);
         
         if (notdone)
            fm=trapz(xs,fs,1);
            xm=ys;
            Hs(i)=trapz(xm,-(fm+1e-14).*log(fm+1e-14));
            notdone=0;
         end
         Mut2ds(i,k)=-H+Hs(i)+Hs(k);
         
         weights(i,k)=Mut2ds(i,k);
         corr=sigma(i,k)/sigma(i,i)^.5/sigma(k,k)^.5;
         
         weights(i,k)=weights(i,k)-lambdaG*.5*log(1-corr^2)+ ...
            +lambdaC*.5*log(1+1e-14-corr^2);
         weights(k,i)=weights(i,k);
         
         
         Mut2ds(k,i)=Mut2ds(i,k);
      end
   end
   
   [etree,weightsum] = MWST(m,weights);
   score=sum(Hs)-log(abs(det(W)));
   % only keeps the first edges
   tree=tree(:,1:param.maxedges);
   for i=1:param.maxedges
      score=score-lambdaT*weights(tree(1,i),tree(2,i));
   end
   
   
   
   
else
   % tree is fixed: a lot less to reevaluate !!!
   adjmat=zeros(m);
   for k=1:size(tree,2)
      adjmat(tree(1,k),tree(2,k))=1;
      adjmat(tree(2,k),tree(1,k))=1;
   end
   toreevaluate=neighbors(adjmat,i);
   notdone=1;
   for k=toreevaluate
      [xs,ys,fs]=kdefft2d(x([k i],:),ngrid,h);
      fs=fs';  
      mu=(fs+1e-14).*log(fs+1e-14);
      H=-trapz(xs,trapz(ys,mu,2),1);
      H2ds(i,k)=H;
      H2ds(k,i)=H2ds(i,k);
      
      if (notdone)
         fm=trapz(xs,fs,1);
         xm=ys;
         Hs(i)=trapz(xm,-(fm+1e-14).*log(fm+1e-14));
         notdone=0;
      end
      Mut2ds(i,k)=-H+Hs(i)+Hs(k);
      Mut2ds(k,i)=Mut2ds(i,k);
   end
   score=sum(Hs)-log(abs(det(W)));
   for k=1:size(tree,2)
      score=score-Mut2ds(tree(1,k),tree(2,k));
   end
   
   score=score+lambdaG*Gaussian_treescore(sigma,tree);
   score=score+lambdaC*correlation_treescore(sigma,tree);
   
   
   
end
