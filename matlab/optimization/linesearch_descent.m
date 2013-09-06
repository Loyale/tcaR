 function [lambda,evalct,fx]= linesearch_descent(fun,x,d,optparam,varargin);

 [lambda,evalct,fx] = linesearch(@directionalfun,optparam,fun,x,d,varargin{:});
