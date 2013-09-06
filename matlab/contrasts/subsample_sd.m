function newperiodo=subsample_sd(periodo,newh);

% only authorize exact subsampling.

[m k oldh]=size(periodo);
newperiodo=zeros(m,m,newh);
for i=1:m
    for j=1:i
        p=reshape(periodo(i,j,:),1,oldh);
        newperiodo(i,j,:)=reshape(subsample(p,newh),1,1,newh);
        if j<i,
            newperiodo(j,i,:)=conj(newperiodo(i,j,:));
        end
        
    end
end

