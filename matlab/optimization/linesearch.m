function [x,evalct,fx] = linesearch(fun,optparam,varargin)
% LINESEARCH - performs line search (algorithm from Numerical Recipes)

scale=optparam.scale_ls;
etafrac=optparam.epsfrac_ls;
etaabs=optparam.epsabs_ls;
Kmax=optparam.maxiter_ls;
% etafrac  = fractional accuracy
% etaabs  = anbsolute accuracy
% scale = first try of bracketing
% Kmax = maximum number of function evaluations for line search
maxdil = 400;                    % Maximum number of dilations allowed
Zeps = etaabs/scale;
eta = etafrac/scale;
if isfield(optparam,'startingstep')
    startingstep=optparam.startingstep;
else
    startingstep=1;
end
%% Phase I: generate a bracketing triple.
%% Triple will have form (0,tm,tr), in golden-mean proportions.
%% Suffixes:  l=left, m=mid, r=right.
%% Assumes that direction vv is a descent direction from point zz.

gold   = (3.0+sqrt(5.0))/2.0;   % Golden ratio expansion factor

tl = 0.0; fl = feval(fun,tl,varargin{:});  % Starting interval has endpoints
%tm = scale; fm = feval(fun,tm,varargin{:});    % tl=0, tm=1=tr.
tm = startingstep; fm = feval(fun,tm,varargin{:});    % tl=0, tm=1=tr.
tr = tm;  fr = fm;

if fr<fl,
    % This is the expected case, where f decreases from t=0 to t=1.
    tr = gold*tr; fr = feval(fun,tr, varargin{:}); ;
    bnum = 1;
    while (fm > min([fl,fr])) & (bnum < maxdil),
        tm = tr;  fm = fr;
        tr = gold*tr; fr = feval(fun,tr, varargin{:}); 
        bnum = bnum+1;
    end;
else
    % This is the strange case, where f increases from t=0 to t=1.
    % Here we shrink the interval instead of expanding it.
    tm = tm/gold; fm = feval(fun,tm,varargin{:}); 
    bnum = 1;
    while (fm > min([fl,fr])) & (bnum < maxdil),
        tr = tm;  fr = fm;
        tm = tm/gold; fm = feval(fun,tm, varargin{:}); ;
        bnum = bnum+1;
    end;
end;
% diagnostic printing
l10tm = log10(tm);
xp = ceil(abs(l10tm))*sign(l10tm);

if fm > min([fl,fr]),
    error('linsrch:  Failed to bracket a minimum');
    evalct = bnum+2;
    lambda = tm;
    return;
end;

%% Phase II:  Apply Brent's Method
%%  This is just typed from Numerical Recipes, modified for Matlab.

% Initializations...
a = tl;  % Left bracket point
b = tr;  % Right bracket point
v = tm;  % Middle bracket point
x = v;   % Point with best fcn value so far
w = v;   % Point with second-best fcn value so far
% (Later, v stores the previous value of w.)
e = 0.0; % Distance moved on the step before last.

t = x; fx = fm;  % Remember fcn value at middle pt of bracket from before.
fv = fx;
fw = fx;


for iter=1:Kmax         % Main program loop.
    
    
    xm = 0.5*(a+b);  % Midpoint of current bracket
    tol1 = eta*abs(x)+Zeps;
    tol2 = 2.*tol1;
    
    % Test for done here.
    if ( abs(x-xm) <= (tol2-0.5*(b-a)) ),
        % Exit with best values.
        lambda = x;
        evalct  = iter+bnum+1;
        return
    end;
    
    if ( abs(e) <= tol1),
        % Step before last was very small, so
        % take a Golden Section Step:
        if (x >= xm),
            e = a-x;
        else
            e = b-x;
        end;
        d = e/gold;
    else
        % Step before last was of decent size, so 
        % construct a trial parabolic fit.
        r = (x-w)*(fx-fv);
        q = (x-v)*(fx-fw);
        p = (x-v)*q - (x-w)*r;
        q = 2.0*(q-r);
        if (q>0.0), p = -p; end;
        q = abs(q);
        etemp = e;
        e = d;
        
        % Test viability of trial fit.
        if (abs(p)>= abs(0.5*q*etemp)) | (p <= q*(a-x)) | ( p>= q*(b-x)) 
            % Parabolic fit is poor, so take golden section step.
            if (x >= xm),
                e = a-x;
            else
                e = b-x;
            end;
            d = e/gold;
        else
            % Parabolic fit is OK, so use it.
            d = p/q;
            u = x+d;
            if (u-a < tol2) | (b-u < tol2), d = tol1*sign(xm-x); end;
        end;
    end;
    
    % Arrive here with increment  d  in hand.
    % Use d to form new x-value, insisting on moving at least tol1.
    if (abs(d) >= tol1),
        u = x+d;
    else
        u = x + tol1*sign(d);
    end;
    
    
    % Evaluate given function at point u
    % (This is the one function evaluation per iteration.)
    fu = feval(fun,u, varargin{:}); 
    if (fu <= fx),
        % New evaluation point  u  is better than best-so-far  x
        if (u >= x), a=x; else b=x; end;   % Contract bracketing interval
        v=w; fv=fw;
        w=x; fw=fx;
        x=u; fx=fu;
    else
        % New evaluation point  u  is worse than best-so-far  x
        if (u < x), a=u; else b=u; end;
        if (fu <= fw) | (w==x),
            v = w; fv = fw;
            w = u; fw = fu;
        elseif(fu<=fv) | (v==x) | (v==w),
            v=u; fv=fu;
        end;
    end;
    
end;    % End of main program loop

evalct = -(iter+bnum+1);
return;
