function [W,tree]=tca_statio(x);
% CLUSTER_TCA - Tree-dependent Component Analysis of stationary time series
%               with default options
%
% for a function with options, seee tca_options.m

n=size(x,2);

fprintf('Tree-dependent Component Analysis - stationary signals + clusters\n')
fprintf('statio...')
[W,tree]=tca_options(x,'contrast','statio','opttype','exhaustive','whitening',1,'display',0,'tolvar',1e-4,'tolval',1e-4);
fprintf('...')
[W,tree]=tca_options(x,'contrast','statio','opttype','localsearch','whitening',1,'W',W,'edgeprior',edgeprior,'display',0,'tolvar',1e-4,'tolval',1e-4);
