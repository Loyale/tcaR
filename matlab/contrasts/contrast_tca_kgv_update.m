function      [score,details]=contrast_tca_kgv_update(x,W,param,i,details);


if (param.fixedtree)
   fixedtree=1;
   tree=param.tree;
else fixedtree=0; end

[m N]=size(x);
x=W'*x;
sigma=x*x'/N;

kkappa=param.kappa;
keta=param.eta;
lambdaG=param.lambdaG;
lambdaC=param.lambdaC;
if (isfield(param,'lambdaT'))
   lambdaT=param.lambdaT;
else
   
   lambdaT=1;
end


Rkappa=details.Rkappa;
Us=details.Us;
Lambdas=details.Lambdas;
Drs=details.Drs;
sizes=details.sizes;
oldstarts=details.starts;
oldsizes=sizes;
weights=details.weights;

%redo the one cholesky decomposition
switch (param.kernel)
      case 'hermite'
     [G,Pvec] =chol_hermite(x(i,:),param.sigma,param.p,N*param.eta); 
   case 'gaussian'
      [G,Pvec] =chol_gauss(x(i,:)/param.sigma,1,N*param.eta); 
         case 'lineargaussian'
      [G,Pvec] =chol_gauss(x(i,:),param.sigma,param.lambdaKG,N*param.eta); 

   case 'poly'
      [G,Pvec] =chol_poly(x(i,:),param.r,param.s,param.d,N*param.eta); 
   end

   
   
   [a,Pvec]=sort(Pvec);
G=centerpartial(G(Pvec,:));

% regularization (see paper for details)
[A,D]=eig(G'*G);
D=diag(D);
indexes=find(D>=N*keta & isreal(D)); %removes small eigenvalues
[newinds,order]=sort(D(indexes));
order=flipud(order);
neig=length(indexes);
indexes=indexes(order(1:neig));  
if (isempty(indexes)), indexes=[1]; end
D=D(indexes);
V=G*(A(:,indexes)*diag(sqrt(1./(D))));
Us{i}=V;
Lambdas{i}=D;
Dr=D;
for j=1:length(D)
   Dr(j)=D(j)/(N*kkappa+D(j));
end
Drs{i}=Dr;
sizes(i)=size(Drs{i},1);

starts=cumsum([1 sizes]);
starts(m+1)=[];

% now creates a new Rkappa, we know that ii is less than jj
newRkappa=eye(sum(sizes));
ii=i;
for i=2:m
   for j=1:i-1
      if ( (j==ii) | (i==ii))
         newbottom=diag(Drs{i})*(Us{i}'*Us{j})*diag(Drs{j});
         newRkappa(starts(i):starts(i)+sizes(i)-1,starts(j):starts(j)+sizes(j)-1)=newbottom;
         newRkappa(starts(j):starts(j)+sizes(j)-1,starts(i):starts(i)+sizes(i)-1)=newbottom';
      else
         newbottom= Rkappa(oldstarts(i):oldstarts(i)+oldsizes(i)-1,oldstarts(j):oldstarts(j)+oldsizes(j)-1);
         newRkappa(starts(i):starts(i)+sizes(i)-1,starts(j):starts(j)+sizes(j)-1)=newbottom;
         newRkappa(starts(j):starts(j)+sizes(j)-1,starts(i):starts(i)+sizes(i)-1)=newbottom';
      end
   end
end

Rkappa=newRkappa;
clear newRkappa;
J=-.5*log(det(Rkappa));
% split, whether we need to calculate the tree or not
if (lambdaT>0)
if (fixedtree)
   for k=1:size(tree,2)
      i=tree(1,k); j=tree(2,k);
      if ((i~=ii) & (j~=ii))
      else
         
         
         Rkappap1=Rkappa(starts(i):starts(i)+sizes(i)-1,starts(j):starts(j)+sizes(j)-1);
         Rkappasm=[eye(sizes(i)) Rkappap1; Rkappap1' eye(sizes(j))];
         weights(i,j)=-.5*log(det(Rkappasm));
         weights(j,i)=weights(i,j);
      end
      
      J=J-lambdaT*weights(i,j);
   end
   J=J+lambdaG*Gaussian_treescore(sigma,tree);
   J=J+lambdaC*correlation_treescore(sigma,tree);
   
   
else
   for i=1:m-1
      for j=i+1:m
         if ((i~=ii) & (j~=ii))
            
         else
            Rkappap1=Rkappa(starts(i):starts(i)+sizes(i)-1,starts(j):starts(j)+sizes(j)-1);
            Rkappasm=[eye(sizes(i)) Rkappap1; Rkappap1' eye(sizes(j))];
            weights(i,j)=-.5*log(det(Rkappasm));
            corr=sigma(i,j)/sigma(i,i)^.5/sigma(j,j)^.5;
            
            weights(i,j)=weights(i,j)-lambdaG*.5*log(1-corr^2)+ ...
               +lambdaC*.5*log(1+1e-14-corr^2);
            
            weights(j,i)=weights(i,j);
         end
         
      end
   end
   details.weights=weights;
   [tree, weightsum] = MWST( m, weights );
      % only keeps the first edges
      tree=tree(:,1:tparam.maxedges);
      for i=1:tparam.maxedges
         J=J-lambdaT*weights(tree(1,i),tree(2,i));
      end

   
   

end
else
   details.weights=[];
end

score=J;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G2=centerpartial(G1)
% CENTERPARTIAL - Center a gram matrix of the form K=G*G'

[N,NG]=size(G1);
G2 = G1 - repmat(mean(G1,1),N,1);






