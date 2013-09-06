% demo script for TCA - 2
% clusters - IID signals
% in case of any bugs or errors, please contact fbach@cs.berkeley.edu

dir=cd;
addpath([dir '/utilities']);
addpath([dir '/contrasts']);
addpath([dir '/optimization']);

load demo_tca2
fprintf('demo script for TCA with clusters - should run in less than 2 minutes\n');
[W,clusters]=cluster_tca(x);
fprintf('estimated clusters (should be two clusters of size two)\n')
for i=1:length(clusters)
clusters{i}
end