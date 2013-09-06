function y=subsample(x,Tsub);
if length(size(x))==2
    T=size(x,2);
    h=T/Tsub;
    filter=[ones(1,h/2) .5 zeros(1,T-h-1) .5 ones(1,h/2-1) ]/h;
    y=ifft(fft(x).*fft(filter));
    y=y(1:h:T);
else
    m1=size(x,1);
    m2=size(x,2);
    T=size(x,3);
    h=T/Tsub;
    filter=[ones(1,h/2) .5 zeros(1,T-h-1) .5 ones(1,h/2-1) ]/h;
    y=x;
    for i=1:m1
        for j=1:m2
            xa=reshape(x(i,j,:),1,T);
            y(i,j,:)=reshape( ifft(fft(xa).*fft(filter)) , 1,1,T);
            
        end
    end
    y=y(:,:,1:h:T);
    
end