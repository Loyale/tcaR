function [W,clusters]=cluster_tca(x);
% CLUSTER_TCA - Tree-dependent Component Analysis of IID signals
%               (with clustering of components)
%               with default options
%
% for a function with options, seee tca_options.m

n=size(x,2);
edgeprior=3*log(n)/n;

fprintf('Tree-dependent Component Analysis - IID signals + clusters\n')
fprintf('ica...')
W=jader(x);
fprintf('cum...')
[W,tree]=tca_options(x,'contrast','cum','opttype','exhaustive','whitening',1,'W',W,'edgeprior',edgeprior,'display',0);
fprintf('...')
[W,tree]=tca_options(x,'contrast','cum','opttype','localsearch','whitening',1,'W',W,'edgeprior',edgeprior,'display',0);
fprintf('kgv...')
[W,tree]=tca_options(x,'contrast','kde','opttype','localsearch','whitening',1,'W',W,'edgeprior',edgeprior,'h',.5,'ngrid',32,'display',0);
fprintf('...')
[W,tree]=tca_options(x,'contrast','kgv','opttype','localsearch','whitening',1,'W',W,'edgeprior',edgeprior,'display',0);
fprintf('kde...')
[W,tree]=tca_options(x,'contrast','kde','opttype','exhaustive','whitening',1,'W',W,'edgeprior',edgeprior,'display',0);
fprintf('...')
[W,tree]=tca_options(x,'contrast','kde','opttype','localsearch','whitening',1,'W',W,'edgeprior',edgeprior,'display',0);

% build clusters
clusters=connected_components(tree2adjmat(tree));