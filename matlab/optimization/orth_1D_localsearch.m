function [score,W]=orth_1D_localsearch(xnorm,W,param,tolx,tolf,nsweeps,detail);
% ORTH_1D_LOCALSEARCH - local search with whiteness constraints

[m n] = size(xnorm);
optparam.scale_ls=pi/8;
optparam.epsfrac_ls=0;
optparam.epsabs_ls=tolx/2;
optparam.maxiter_ls=40;

[score,details]=contrast_tca(xnorm,W,param);
for k=1:nsweeps
    errorx=0;
    errorf=0;
    oldscore=score;
    for i=2:m
        for j=1:i-1
            oldtree=details.tree;
            param.fixedtree=1;
            param.tree=details.tree;
            [theta,evalct,fx] = linesearch(@orth_tca,optparam,xnorm,W,param,details,i,j);
            c = [cos(theta) sin(theta); ...
                    -sin(theta) cos(theta);   ];
            W([i j],:)= c * W([i j],:);
            param.fixedtree=0;
            [score,details]=contrast_tca_multupdate(xnorm,W,param,[i j],details);
           % [score,details]=contrast_tca(xnorm,W,param);
            newtree=details.tree;
          %  fprintf('i=%d j=%d itheta1=%f oldscore=%f newscore=%f newtree=%d\n',i,j,theta,oldscore,score,1-compare_tree(m,newtree,oldtree));
                errorx=errorx+abs(theta)^2;
        end
    end
      if detail,  fprintf('sweep %d, oldscore=%f newscore=%f deltax=%f\n',k,oldscore,score,sqrt(errorx)); end
    errorf=errorf+abs(score-oldscore);
    if sqrt(errorx)<tolx | errorf<tolf, return; end
end