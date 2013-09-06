% demo script for TCA - 3
% clusters - IID signals
% in case of any bugs or errors, please contact fbach@cs.berkeley.edu

dir=cd;
addpath([dir '/utilities']);
addpath([dir '/contrasts']);
addpath([dir '/optimization']);


load demo_tca3
fprintf('demo script for TCA with clusters and stationary signals - should run in less than 1 minute\n');
[W,clusters]=cluster_tca_statio(x);
fprintf('estimated clusters (should be two clusters of size two)\n')
for i=1:length(clusters)
clusters{i}
end