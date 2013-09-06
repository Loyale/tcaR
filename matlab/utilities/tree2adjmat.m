function adjmat=tree2adjmat(edges,m);
if nargin<2,m=size(edges,2)+1; end
adjmat=zeros(m);
for i=1:size(edges,2)
   adjmat(edges(1,i),edges(2,i))=1;
   adjmat(edges(2,i),edges(1,i))=1;
end
