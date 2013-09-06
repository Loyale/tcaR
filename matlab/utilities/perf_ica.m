function [perf,assignement]=perf_ica(W,A);



m=size(W,1);

P=abs(W*A)./repmat(sum(abs(W*A),2),1,m);
assignement=zeros(1,m);
perf=0;
for i=1:m
    [a,b]=max(P(i,:));
    perf=perf+sum(P(i,:))/a;
    assignement(i)=b;
end

perf=(perf-m)/(m-1)/m;


return

P=W*A;
% build assignement cost matrix
j=1;
C=zeros(m);
for i=1:m;
    for k=1:m
    C(i,k)=sum(abs(P(:,k)))/abs(P(i,k))-1;
end
end
C
[assignement,perf]=hungarian(C);
