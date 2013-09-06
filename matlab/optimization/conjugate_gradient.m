function [x,fx, nfun, ngfun] = conjugate_gradient(fun,gfun,x,optparam,varargin)
% CONJUGATE_GRADIENT : performs conjugate gradient optimization

if ~isfield(optparam,'global_cg'), optparam.global_cg=0; end
eps_x=optparam.epsx_cg;
eps_f=optparam.epsf_cg;
Kmax=optparam.maxiter_cg;
display=optparam.display;
optparam.epsfrac_ls=eps_x;
optparam.epsabs_ls=eps_x;
empgradtype=optparam.empgradtype;


if ~isa(gfun,'function_handle')
    % use empirical gradients
    empgrad=1;
else
    empgrad=0;
end
xtyp = ones(size(x));
ftyp = sqrt(eps);

% Initialization.
xo = x;
k  = 0;                       
dx = 10*eps_x;
df = 10*eps_f;
oldt=1;
nfun = 0;
ngfun = 0;
while or(dx>eps_x,df>eps_f*16) &  or(df>eps_f,dx>eps_x*16)  & (k<Kmax) & (dx>eps_x/8) & (df>eps_f/8)
    [fx, parameters] = feval(fun,x,'global',[],varargin{:});  % Function value at initial point.
    nfun = nfun + 1;
    if display,  fprintf('g...'); end
    if empgrad
        % compute empirical gradient
        scaling = 1;
        v = compute_empgrad(fun,x,fx,scaling,empgradtype,'local',parameters,varargin{:});
    else
        v  = feval(gfun,x,varargin{:}); % Gradient at initial point (row vector).
    end
    ngfun = ngfun + 1;
    if display,  fprintf('done. '); end
    
    % Descent direction
    if k==0,
        d = -v; 
        b=0;
    else
        d = -v + ((v-vo)'*v) / (vo'*vo) * d;               % Polak-Ribiere update rule
        if b'*v >= 0, % make sure this is a descent direction
            d = -v;
        end
    end;
    if display,  fprintf('b=%f - ',b); end
    
    
    % Do line search from base point x in direction d.
    optparam.scale_ls=norm(d)+eps_f;
    % record the previous amount of line search and start with it??
    optparam.startingstep=oldt;
    
    optparam.maxiter_ls=100;
    fxo=fx;
    if display,  fprintf('fx=%f - gnorm=%f - ',fx,norm(d)); end
    if ~optparam.global_cg,
        [t,evals,fx] = linesearch_descent(fun,x,d,optparam,'local',parameters,varargin{:});
    else
        [t,evals,fx] = linesearch_descent(fun,x,d,optparam,'global',[],varargin{:});
    end
    nfun = nfun + evals;
    if display, fprintf('nev=%d - df=%f - dx=%f\n',evals,fxo-fx,norm(t*d)); end
    xo = x;         
    vo = v;         
    x = x + t*d; 
    oldt=t;
    dx = norm(x-xo);
    df = norm(fx-fxo);
    k = k+1;
    
end;      

evalct = k;
return
