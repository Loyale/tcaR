function [W,tree,details]=tca_options(x,varargin)

% TCA_OPTIONS - runs the tca algorithm on the m x n data x
%               m = number of dimensions
%               n = number of samples
%
% OPTIONS:
% 'whitening'      : 1 -> whitening, 0 -> no whitening
% 'opttype'        : 'exhaustive' or 'localsearch'
% 'tolvar'         : precision in variable W
% 'tolval'         : precision in value
% 'lambdaG'        : weight of Gaussian contrast (default = 0)
% 'lambdaT'        : weight of mutual informatin terms (default = 1)
% 'lambdaC'        : weight of correlation penalization term (default = 0.05)
% 'maxedges'       : maximum number of edges allowed in the tree
% 'edgeprior'      : penanlty term for each edge (non zero to find clusters)
% 'nsweeps'        : maximum number of sweeps
% 'searchgrid'     : number of points in the grid for exhaustive search
% 'W'              : initial guess
%
% 'contraste'      : 'cum' 'kgv' 'kde' 'statio'
% 'h' 'ngrid'      : parameter for the contrast 'kde'
% 'kappa'          : refularization parameter for contrast 'kgv'
% 'kernel'         : 'gaussian' 'hermite' 'poly'
% 'r' 's'          : parameters for kernel 'poly'
% 'sigma' 'd'      : parameters for kernel 'hermite'
% 'sigma'          : parameter for kernel 'gaussian'
% 'display'        : 1 -> show information, 0-> silent

% default parameters
display =1;
contraste = 'cum';
[m n] = size(x);
kernel='gaussian';
p=4;
sigma=1;
kappa=.01;
lambdaG=0;
lambdaT=1;
lambdaC=0.05;
h=.25;
ngrid=64;
tolvar=4e-4;
tolval=4e-4;
maxedges=m-1;
whitening=0;
opttype = 'localsearch';
edgeprior = 0;
W=[];
r=1;
s=1;
d=4;
Wgold=[];
searchgrid=16;
nsweeps=m*m;

% reading parameters
args = varargin;
nargs = length(args);
for i=1:2:nargs
    switch args{i},
        case 'edgeprior',  edgeprior = args{i+1};
        case 'searchgrid',  searchgrid = args{i+1};
        case 'tolvar',      tolvar = args{i+1};
        case 'tolval',      tolval = args{i+1};
        case 'opttype',     opttype = args{i+1};
        case 'whitening',   whitening = args{i+1};
        case 'contrast',    contraste = args{i+1};
        case 'h',           h = args{i+1};
        case 'r',           r = args{i+1};
        case 's',           s = args{i+1};
        case 'd',           d = args{i+1};
        case 'kappa',       kappa = args{i+1};
        case 'sigma',       sigma = args{i+1};
        case 'ngrid',       ngrid = args{i+1};
        case 'lambdaG',     lambdaG = args{i+1};
        case 'lambdaC',     lambdaC = args{i+1};
        case 'lambdaT',     lambdaT = args{i+1};
        case 'maxedges',    maxedges = args{i+1};
        case 'kernel',      kernel = args{i+1};
        case 'W',           W = args{i+1};
        case 'Wgold',       Wgold = args{i+1};
        case 'nsweeps',     nsweeps = args{i+1};
        case 'display',     display = args{i+1};
    end
end

param.contraste=contraste;
param.kernel=kernel;
param.p=p;
param.r=r;
param.s=s;
param.d=d;
param.sigma=sigma;
param.kappa=kappa;
param.eta=kappa*1e-2;
param.lambdaG=0;
param.lambdaT=lambdaT;
param.lambdaC=lambdaC;
param.h=h;
param.ngrid=ngrid;
param.maxedges=maxedges;
param.edgeprior=edgeprior;


% first compute the mean and covariance
mu = mean(x,2);
sigma = 1 / n * ( x - repmat(mu,1,n) ) * ( x - repmat(mu,1,n) )';
G = chol( sigma );
Wc = inv(G');
xnorm = Wc * ( x - repmat(mu,1,n) ) ;

if strcmp(param.contraste,'statio')
    % for the statio contraste, compute the joint spectral density
    
    % compute optimal marginal smoothing widths
    for i=1:m
        [sdloc,dfloc,rloc]=spectral_density(xnorm(i,:),'display',0,'mstype','aic');
        dflocs(i)=dfloc;
        rlocs(i)=rloc;
    end
    rmean=mean(rlocs);
    df=mean(dflocs);
    
    % smooth periodogram
    hsd=min([ 2^( round(log2(df*4))), n]);
    [Wfft,df]=generate_smoothingwin(n,'gaussian',rmean);
    smoothingprior=1e-3;
    periodo=myperiodogram(xnorm);
    for t=1:n
        periodo(:,:,t)=periodo(:,:,t)+smoothingprior*eye(m);
    end
    sd=zeros(m,m,n);
    for i=1:m
        for j=1:m
            peri=reshape(periodo(i,j,:),[ 1 n]);
            peri=ifft(Wfft.*fft(peri))/2/pi;
            sd(i,j,:)=reshape(peri,[ 1 1 n ]);
            
        end
    end
    sd=subsample_sd(sd,hsd);
    param.periodo=sd;
    param.df=df;
end

% initialization
if isempty(W), 
    W = rand_orth(m); 
else 
    W = W * inv( Wc );
    if whitening,
        [U,S,V]=svd(W); W=U*V'; 
    end
end



if whitening
    switch opttype
        case 'exhaustive'
           if display, fprintf('SWEEP 1\n'); end
            [score,W]=orth_1D_exhaustive(xnorm,W,param,round(searchgrid/2),nsweeps,display);
           if display, fprintf('SWEEP 2\n'); end
            [score,W]=orth_1D_exhaustive(xnorm,W,param,round(searchgrid),nsweeps,display);
        case 'localsearch'
            [score,W]=orth_1D_localsearch(xnorm,W,param,tolvar,tolval,nsweeps,display);
    end
    
else
    switch opttype
        case 'exhaustive'
            [score,W]=spheres_2D_exhaustive(xnorm,W,param,round(searchgrid),nsweeps,display);
        case 'localsearch'
            [score,W]=spheres_2D_localsearch(xnorm,W,param,tolvar,tolval,nsweeps,display);
    end
end


% computes the full tree (without any penalty)
param.maxedges=m-1;
[scoregold,details]=contrast_tca(xnorm,W,param);
treefull=details.tree;
param.maxedges=maxedges;
[scoregold,details]=contrast_tca(xnorm,W,param);
tree=details.tree;
details.score=scoregold;
details.treefull=treefull;
W = W * Wc;

% make sure that the 'leaf mixing' invariance is normalized
W=normalize_tcaresult(x,W,tree);
