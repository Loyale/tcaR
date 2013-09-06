% demo script for TCA - 1
% IID signals
% in case of any bugs or errors, please contact fbach@cs.berkeley.edu

dir=cd;
addpath([dir '/utilities']);
addpath([dir '/contrasts']);
addpath([dir '/optimization']);
load demo_tca1
fprintf('demo script for TCA - should run in less than 2 minutes\n');
[W,treelearn]=tca(x);
[perf,assignement]=perf_ica(W,A);
fprintf('performance in W = %f \n', perf);
fprintf('original tree\n')
disp((tree));
fprintf('estimated tree\n')
disp(assignement(treelearn));