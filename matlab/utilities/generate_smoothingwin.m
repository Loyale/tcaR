function [Wfft,DF]=generate_smoothingwin(T,smoothingtype,rsmooth);

js=[0:T/2 -T/2+1:-1];
M=(T/2+1);

switch smoothingtype,
    case 'daniell'
        Wfft=js;
        Wfft(2:end)=sin(js(2:end)*pi/rsmooth)./(js(2:end)*pi/rsmooth);
        Wfft(1)=1;
        % pdaniell=indices(smoothingwidth);
        %  W=[ones(1,pdaniell+1) zeros(1,T-2*pdaniell-1) ones(1,pdaniell)]/(2*pdaniell+1);
        %  DF=W(1)*M;
        
        %Wfft=real(fft(W));
        DF=sum(Wfft)/T*M;
                DF=sum(Wfft);

    case 'bartlett'
        for t=1:T
            if abs(js(t)/rsmooth) >1 , Wfft(t)=0;
            else  Wfft(t)=1-abs(js(t)/rsmooth); end
        end
        DF=sum(Wfft)/T*M;
        DF=sum(Wfft);
    case 'gaussian'
        Wfft=exp(-.5*(js/rsmooth).^2); 
        DF=sum(Wfft)/T*M;
                DF=sum(Wfft);

    case 'spline2'
        Wfft=1./(1+(js/rsmooth).^2); 
        DF=sum(Wfft)/T*M;
                DF=sum(Wfft);

    case 'spline4'
        Wfft=1./(1+(js/rsmooth).^4); 
        DF=sum(Wfft)/T*M;
                DF=sum(Wfft);

        
end

