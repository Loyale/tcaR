function [score,W]=spheres_2D_localsearch(xnorm,W,param,tolx,tolf,nsweeps,detail);
[m n] = size(xnorm);
optparam.epsx_cg=tolx/2;
optparam.epsf_cg=tolf/2;
optparam.maxiter_cg=80;
optparam.display=0;
optparam.empgradtype=1;
param.fixedtree=0;
[score,details]=contrast_tca(xnorm,W,param);
finished=zeros(m);
for k=1:nsweeps
    oldscore=score;
    errorx=0;
    errorf=0;
    
    for i=2:m
        for j=1:i-1
            oldtree=details.tree;
            param.fixedtree=1;
            param.tree=details.tree;
            B = W([i j],:) * W([i j],:)' ;
            assert(abs(B(1,1)-1)<1e-8 & abs(B(2,2)-1)<1e-8,'row of W is not unit norm');
            D = chol(B);
            c = [1 0] * D';
            theta10 = angle(c(1)+sqrt(-1)*c(2));
            c = [0 1] * D' ;
            invD = inv( D' );
            theta20 = angle(c(1)+sqrt(-1)*c(2));
            
            [theta,fx]= conjugate_gradient(@sphere_tca,0,[0;0],optparam,xnorm,W,param,details,i,j,invD,theta10,theta20);
            theta1=theta(1);
            theta2=theta(2);
            
            
            c = [cos(theta1+theta10) sin(theta1+theta10); ...
                    cos(theta2+theta20) sin(theta2+theta20)];
            a = c * invD;
            W([i j],:)=a * W([i j],:);
            %            [score,details]=contrast_tca_2update(xnorm,W,param,i,j,details);   
            param.fixedtree=0;
            [score,details]=contrast_tca_multupdate(xnorm,W,param,[i j],details);
            
            newtree=details.tree;
            
            
        %    fprintf('i=%d j=%d itheta1=%f itheta2=%f oldscore=%f newscore=%f newtree=%d\n',i,j,theta1,theta2,oldscore,score,1-compare_tree(m,newtree,oldtree));
        end
    end
    errorx=errorx+.5*abs(theta1)^2+.5*abs(theta2)^2;
    errorf=errorf+abs(score-oldscore);
     if detail,   fprintf('sweep %d, oldscore=%f newscore=%f\n',k,oldscore,score); end
    if sqrt(errorx)<tolx & errorf<tolf, return; end
end