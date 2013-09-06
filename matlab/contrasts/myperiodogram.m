function periodo=periodogram(x);

% Compute the periodogram of a multivariate time series
[m N]=size(x);
means=mean(x,2);
x=x-repmat(means,1,N);


xfft=x;
for i=1:m
    xfft(i,:)=fft(x(i,:));
end
xfft=xfft/sqrt(N);


periodo=zeros(m,m,N);

for i=1:m
    for j=1:m
        temp=xfft(i,:).*conj(xfft(j,:));
        periodo(i,j,:)=reshape(temp(1:N),1,1,N);
    end
end
