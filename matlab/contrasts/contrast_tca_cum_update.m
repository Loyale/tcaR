function      [score,details]=contrast_tca_cum_update(x,W,param,i,details);

if (param.fixedtree)
   fixedtree=1;
   tree=param.tree;
else fixedtree=0; end

[m N]=size(x);
x=W'*x;
sigma=x*x'/N;
Mut2ds=details.Mut2ds;
Hs=details.Hs;
lambdaC=param.lambdaC;
lambdaG=param.lambdaG;


cum4=zeros(1,m);
cum4=sum(x.^4,2)/N;
cum3=sum(x.^3,2)/N;
Hs=-1/48 * (cum4 - 3).^2 - 1/12 * cum3.^2; % + 5/8 * cum3.^2 .* cum4 + 1/16 * (cum4 - 3).^3;
%Hs=zeros(1,m);

if (~fixedtree)
   notdone=zeros(1,m);
   notdone=1;
   for k=1:m
      if (k~=i)
         
                  y=x([i k],:);
         sigy=sigma([i k],[i k]);
         y=inv(sqrtm(sigy))*y;
         Mut2ds(i,k)=-.5*log(1-sigma(i,k)^2) ...
            + 1/12*( sum( (y(1,:).^3).*y(2,:),2)/N )^2 ...
            + 1/8*(  sum( (y(1,:).^2).*(y(2,:).^2),2)/N -1)^2 ...
            + 1/12*( sum( (y(2,:).^3).*y(1,:),2)/N)^2 ...
            + 1/48*( (sum(y(1,:).^4)/N-3)^2 + (sum(y(2,:).^4)/N-3)^2) ...
                        + 1/4*( sum( (y(2,:).^2).*y(1,:),2)/N)^2 ...
            + 1/4*( sum( (y(1,:).^2).*y(2,:),2)/N)^2 ...
            + 1/12*( (sum(y(1,:).^3)/N)^2 + (sum(y(2,:).^3)/N)^2) ...
            + Hs(i)+Hs(k);


      %   Mut2ds(i,k)=-.5*log(1-sigma(i,k)^2);
      %   1 ...
      %      + 1/12*( sum( (x(i,:).^3).*x(k,:),2)/N-3*sigma(i,k))^2 ...
      %      + 1/8*(  sum( (x(i,:).^2).*(x(k,:).^2),2)/N -1 -2*sigma(i,k)^2)^2 ...
      %      + 1/12*( sum( (x(k,:).^3).*x(i,:),2)/N-3*sigma(i,k))^2;
         
         
         weights(i,k)=Mut2ds(i,k);
         corr=sigma(i,k)/sigma(i,i)^.5/sigma(k,k)^.5;
         weights(i,k)=weights(i,k)-lambdaG*.5*log(1-corr^2)+ ...
            +lambdaC*.5*log(1+1e-14-corr^2);
         weights(k,i)=weights(i,k);
         
         
         Mut2ds(k,i)=Mut2ds(i,k);
      end
   end
   
   [etree,weightsum] = MWST(m,weights);
   score=sum(Hs)-(1+lambdaG)*log(abs(det(W)));
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
      
                        y=x([i k],:);
         sigy=sigma([i k],[i k]);
         y=inv(sqrtm(sigy))*y;
         Mut2ds(i,k)=-.5*log(1-sigma(i,k)^2) ...
            + 1/12*( sum( (y(1,:).^3).*y(2,:),2)/N )^2 ...
            + 1/8*(  sum( (y(1,:).^2).*(y(2,:).^2),2)/N -1)^2 ...
            + 1/12*( sum( (y(2,:).^3).*y(1,:),2)/N)^2 ...
            + 1/48*( (sum(y(1,:).^4)/N-3)^2 + (sum(y(2,:).^4)/N-3)^2) ...
            + 1/4*( sum( (y(2,:).^2).*y(1,:),2)/N)^2 ...
            + 1/4*( sum( (y(1,:).^2).*y(2,:),2)/N)^2 ...
            + 1/12*( (sum(y(1,:).^3)/N)^2 + (sum(y(2,:).^3)/N)^2) ...
            + Hs(i)+Hs(k);


    %  Mut2ds(i,k)=-.5*log(1-sigma(i,k)^2);
    %  1 ...
    %        + 1/12*( sum( (x(i,:).^3).*x(k,:),2)/N-3*sigma(i,k))^2 ...
    %        + 1/8*(  sum( (x(i,:).^2).*(x(k,:).^2),2)/N -1 -2*sigma(i,k)^2)^2 ...
    % + 1/12*( sum( (x(k,:).^3).*x(i,:),2)/N-3*sigma(i,k))^2;
      Mut2ds(k,i)=Mut2ds(i,k);
   end
   score=sum(Hs)-log(abs(det(W)));
   for k=1:size(tree,2)
      score=score-Mut2ds(tree(1,k),tree(2,k));
   end
   
   score=score+lambdaG*Gaussian_treescore(sigma,tree);
   score=score+lambdaC*correlation_treescore(sigma,tree);
   
   
   
end
