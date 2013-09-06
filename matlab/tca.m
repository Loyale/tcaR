function [W,tree]=tca(x);
% TCA - Tree-dependent Component Analysis of IID signals
%       with default options
%
% for a function with options, seee tca_options.m


fprintf('Tree-dependent Component Analysis - IID signals\n')
fprintf('ica...')
W=jader(x);
fprintf('cum...')
[W,tree]=tca_options(x,'contrast','cum','opttype','exhaustive','whitening',1,'W',W,'display',0);
fprintf('...')
[W,tree]=tca_options(x,'contrast','cum','opttype','localsearch','whitening',1,'W',W,'display',0);
fprintf('kgv...')
[W,tree]=tca_options(x,'contrast','kde','opttype','localsearch','whitening',0,'W',W,'h',.5,'ngrid',32,'display',0);
fprintf('...')
[W,tree]=tca_options(x,'contrast','kgv','opttype','localsearch','whitening',0,'W',W,'display',0);
fprintf('kde...')
[W,tree]=tca_options(x,'contrast','kde','opttype','localsearch','whitening',0,'W',W,'display',0);
