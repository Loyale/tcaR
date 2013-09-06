function    v = compute_empgrad(fun,x,fx,scaling,empgradtype,varargin);

% COMPUTE_EMPGRAD - computes the empirical gradient
% empgradtype  : 1 - d function evaulations
%                2 - 2*d functions evaulations

d = size(x,1);
dr = 1e-4;
v = zeros(d,1);
y0 = x;
fy0 = fx;
for i=1:d
    y1 = y0; 
    y1(i) = y1(i) + dr * scaling;
    fy1 = feval(fun,y1,varargin{:}); 
    if empgradtype==2
        y2 = y0;
        y2(i) = y2(i) - dr * scaling;
        fy2 = feval(fun,y2,varargin{:}); 
        v(i) = ( fy1 - fy2 ) / 2 / scaling / dr ;
    else
       v(i) = ( fy1 - fy0 ) / scaling / dr ;
    end
end
