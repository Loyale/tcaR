function f=directionalfun(t,fun,x,d,varargin);
f=feval(fun,x+t*d,varargin{:});
