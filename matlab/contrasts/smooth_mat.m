function [periodo,DF]=smooth_mat(periodo,smoothingtype,rsmooth);

m=size(periodo,1);
T=size(periodo,3);
[Wfft,DF]=generate_smoothingwin(T,smoothingtype,rsmooth);


% smoothes
for i=1:m
    for j=1:i
        %   smoothedperiodo(i,j,:)=reshape(ifft(fft(W).*fft(squeeze(periodo(i,j,:)).')),1,1,T);  
        periodo(i,j,:)=reshape(ifft(Wfft.*fft(squeeze(periodo(i,j,:)).')),1,1,T)/2/pi;  
        if j<i
        periodo(j,i,:)=conj(periodo(i,j,:));
        end
    end
end
