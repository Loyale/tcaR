function W=normalize_tcaresult(x,W,tree)
% normalizing W
% we only assume that the components have unit variance
m=size(x,1);
n=size(x,2);
mu = mean(x,2);
sigma = 1 / n * ( x - repmat(mu,1,n) ) * ( x - repmat(mu,1,n) )';
adjmat=zeros(m);
for k=1:size(tree,2)
    adjmat(tree(1,k),tree(2,k))=1;
    adjmat(tree(2,k),tree(1,k))=1;
end
for leaf=1:m
    if (length(neighbors(adjmat,leaf))==1)
         % for each leaf
         parent=neighbors(adjmat,leaf);
         sigmaloc=W([leaf parent],:)*sigma*W([leaf parent],:)';
         assert( abs(sigmaloc(1,1)-1)<1e-6 &  abs(sigmaloc(2,2)-1)<1e-6 , 'not unit variance' );
         rho=sigmaloc(1,2);
        W(leaf,:)=1/sqrt(1-rho^2)*W(leaf,:)- rho/sqrt(1-rho^2)*W(parent,:);

    end
end
