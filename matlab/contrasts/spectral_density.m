function [periodo,DF,rsmooth]=spectral_density(samples,varargin);

[m T]=size(samples);

args = varargin;
nargs = length(args);
smoothingtype='gaussian';
smoothingwidth=10;
smoothingprior=1e-3;
mstype='aic';
ms=1;
% regparam is a mix of aic and bic, for small T, aic, for large T, bic
alpha=max([ 1 log(T)/log(256)]);
regparam=(.5*log(T))^( .5 );
aa=max([2*ceil(T^.3), 2*16]);
hh=ceil(sqrt(2*pi)*T^.3);
searchgrid=1;
bb=floor(2*T^.8);
per0=[];
band=[];
display=1;
for i=1:2:nargs
    switch args{i},
        case 'band',   band = args{i+1};
        case 'per0',   per0 = args{i+1};
        case 'a',   aa = args{i+1};
        case 'b',   bb = args{i+1};
        case 'h',   hh = args{i+1};
        case 'searchgrid',   searchgrid = args{i+1};
        case 'smoothingtype',   smoothingtype = args{i+1};
        case 'bicaic',          bicaic = args{i+1};
        case 'smoothingprior',       smoothingprior = args{i+1};
        case 'mstype',       mstype = args{i+1};
        case 'ms',       ms = args{i+1};
        case 'regparam',       regparam = args{i+1};
        case 'display',       display = args{i+1};
    end
end
hh=hh*searchgrid;
indices=aa:hh:bb;


periodo=myperiodogram(samples);
for t=1:T
    periodo(:,:,t)=periodo(:,:,t)+smoothingprior*eye(m);
end
if ~isempty(per0),   periodo(:,:,1)=per0; end;

smoothedperiodo=periodo;
logl=[];
df=[];
M=(T/2+1);

if ms,
    if display,
        fprintf('bandwidth selection for smoothing periodogram\n');
        fprintf('number of steps=%d, from %d to %d', round((bb-aa)/hh),aa,bb);
    end
    % impose minimal smoothing
    for smoothingwidth=1:length(indices)
        
        rsmooth=T/indices(smoothingwidth);
        js=[0:T/2 -T/2+1:-1];
        if display
            if mod(smoothingwidth,50)==1, fprintf('\n%d',smoothingwidth-1); end
            fprintf('.');
        end
        % builds smooothing window
        [smoothedperiodo,DF]=smooth_mat(periodo,smoothingtype,rsmooth);
        
        % computes the deviance
        deviance=0;
        if m==1
            smperiodo=squeeze(smoothedperiodo);
            perio=squeeze(periodo);
            if isempty(band)
                deviance=sum(log(smperiodo(2:T/2))+perio(2:T/2)/2/pi./smperiodo(2:T/2));
                deviance=deviance+.5*log(smperiodo(1))+.5*perio(1)/2/pi/smperiodo(1);
                deviance=deviance+.5*log(smperiodo(T/2+1))+.5*perio(T/2+1)/2/pi/smperiodo(T/2+1);
                logl(smoothingwidth)=real(deviance)/T;
                df(smoothingwidth)=DF;
                
            else
                deviance=sum(log(smperiodo(band))+perio(band)/2/pi./smperiodo(band));
                
                logl(smoothingwidth)=real(deviance)/length(band)/2;
                df(smoothingwidth)=DF;
                
                
            end
            
        else
            if isempty(band)
                for t=2:(T/2+1)
                    deviance=deviance+log(det(smoothedperiodo(:,:,t)))+trace(periodo(:,:,t)/2/pi/smoothedperiodo(:,:,t));
                end
                deviance=deviance+.5*log(det(smoothedperiodo(:,:,1)))+.5*trace(periodo(:,:,1)/2/pi/smoothedperiodo(:,:,1));
                deviance=deviance+.5*log(det(smoothedperiodo(:,:,T/2+1)))+.5*trace(periodo(:,:,T/2+1)/2/pi/smoothedperiodo(:,:,T/2+1));
                logl(smoothingwidth)=real(deviance)/T;
                df(smoothingwidth)=DF;
                
            else
                for t=1:length(band)
                    deviance=deviance+log(det(smoothedperiodo(:,:,band(t))))+ ...
                        trace(periodo(:,:,band(t))/2/pi/smoothedperiodo(:,:,band(t)));
                end
                logl(smoothingwidth)=real(deviance)/length(band)/2;
                df(smoothingwidth)=DF;
                
                
            end
        end
        
    end
        clear smoothperiodo;

    % [abic,bbic]=min(df/2/M*log(M)*m*m+logl);
    [abic,bbic]=min(df/T*2*log(T)*m*m/2+logl);
    % [aaic,baic]=min(df/M*m*m+logl);
    [aaic,baic]=min(df/T*m*m/2+logl);
    %  [aaicbic,baicbic]=min(df/2/M*(log(M)-log(M)/M^.125 )*m*m+logl);
    % [aaicbic,baicbic]=min(df/2/T*(log(T)-log(T)/T^.125 )*m*m+logl);
    [aaicbic,baicbic]=min(df/T*m*m/2*regparam+logl);
    % [abic,bbic]=min(df/2/M*log(M)*m*(m+1)/2+logl);
    % [aaic,baic]=min(df/M*m*(m+1)/2+logl);
    % [aaicbic,baicbic]=min(df/2/M*(log(M)-log(M)/M^.125 )*m*(m+1)/2+logl);
        switch mstype,
            case 'aic', smoothingwidth=baic; score=aaic;   if display,  plot(df/T*m*m/2+logl); end
            case 'bic', smoothingwidth=bbic; score=abic;  if display,  plot(df/T/2*log(T)*m*m/2+logl); end
            case 'bicaic', smoothingwidth=baicbic; score=aaicbic;  if display,  plot(df/T*m*m/2*regparam+logl); end
        end
        
            if display,

        pause(.1);
        
        %   plot(df/2/M*log(M)*m*(m+1)/2+logl)
        fprintf('\n finished, sw=%d, bw=%d, score=%f\n',smoothingwidth,indices(smoothingwidth),score);
    end
    bw=indices(smoothingwidth);
    rsmooth=T/indices(smoothingwidth);
    clear smoothperiodo;
    [periodo,DF]=smooth_mat(periodo,smoothingtype,rsmooth);
    
end

