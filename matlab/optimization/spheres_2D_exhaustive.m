function [score,W]=spheres_2D_exhaustive(xnorm,W,param,searchgrid,nsweeps,detail);
[m n] = size(xnorm);

[score,details]=contrast_tca(xnorm,W,param);
finished=zeros(m);
for k=1:nsweeps
    for i=2:m
        for j=1:i-1
                        oldtree=details.tree;
            param.fixedtree=1;
            param.tree=details.tree;
            oldscore=score;
            % only optimizes with respect to i and j, while keeping fixed all
            % other ones
            theta1 = (0:searchgrid-1)*pi/searchgrid;
            theta2 = (0:searchgrid-1)*pi/searchgrid;
            B = W([i j],:) * W([i j],:)' ;
            assert(abs(B(1,1)-1)<1e-8 & abs(B(2,2)-1)<1e-8,'row of W is not unit norm');
            D = chol(B);
            c = [1 0] * D';
            theta10 = angle(c(1)+sqrt(-1)*c(2));
            c = [0 1] * D' ;
            invD = inv( D' );
            theta20 = angle(c(1)+sqrt(-1)*c(2));
            for itheta1=1:length(theta1)
                for itheta2=1:length(theta2)
                    W1=W;
                    % builds new W, such that for theta1=theta2=0, we don't
                    % move
                    c = [cos(theta1(itheta1)+theta10) sin(theta1(itheta1)+theta10); ...
                            cos(theta2(itheta2)+theta20) sin(theta2(itheta2)+theta20)];
                    a = c * invD;
                    W1([i j],:)=a * W1([i j],:);
                    %  W1([i j],:)
                    if abs(norm(W1(i,:)-W1(j,:)))<1e-8 || abs(norm(W1(i,:)+W1(j,:)))<1e-8,
                        score = Inf;
                    else
                        [score,details]=contrast_tca_multupdate(xnorm,W1,param,[i,j],details);
                    end
                    scoreloc(itheta1,itheta2) = score;
                end
            end
            % finds minimum and update W
            [a,b]=min(scoreloc(:));
            c=ind2subv([length(theta1) length(theta2)], b);
            itheta1=c(1); itheta2=c(2);
            c = [cos(theta1(itheta1)+theta10) sin(theta1(itheta1)+theta10); ...
                    cos(theta2(itheta2)+theta20) sin(theta2(itheta2)+theta20)];
            a = c * invD;
            W([i j],:)=a * W([i j],:);
            param.fixedtree=0;

            [score,details]=contrast_tca(xnorm,W,param);
 if detail,           fprintf('i=%d j=%d itheta1=%d itheta2=%d oldscore=%f newscore=%f \n',i,j,itheta1-1,itheta2-1,oldscore,score); end
            if abs(oldscore-score)<1e-6
                finished(i,j)=1;
                if sum(finished(:))==m*(m-1)/2, return; end
            else
                finished=zeros(m);
            end
        end
    end
end