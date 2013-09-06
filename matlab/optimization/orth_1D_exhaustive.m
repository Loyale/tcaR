function [score,W]=orth_1D_exhaustive(xnorm,W,param,searchgrid,nsweeps,detail);
% ORTH_1D_EXHAUSTIVE - exhaustive search with whiteness constraints

[m n] = size(xnorm);

[score,details]=contrast_tca(xnorm,W,param);
finished=zeros(m);
for k=1:nsweeps
    reallyoldscore=score;
    for i=2:m
        for j=1:i-1
            oldtree=details.tree;
            
            oldscore=score;
            param.fixedtree=0;
            param.tree=details.tree;
            
            % only optimizes with respect to i and j, while keeping fixed all
            % other ones
            theta1 = (0:searchgrid-1)*pi/searchgrid/2;
            
            B = W([i j],:) * W([i j],:)' ;
            assert(norm(B-eye(2))<1e-6,'W is not orthogonal');
            for itheta1=1:length(theta1)
                W1=W;
                c = [cos(theta1(itheta1)) sin(theta1(itheta1)); ...
                        -sin(theta1(itheta1)) cos(theta1(itheta1));   ];
                W1([i j],:)= c * W1([i j],:);
                
                [score,details]=contrast_tca_multupdate(xnorm,W1,param,[i j],details);
                scoreloc(itheta1) = score;
            end
            % finds minimum and update W
            [a,b]=min(scoreloc);
            itheta1=b; 
            c = [cos(theta1(itheta1)) sin(theta1(itheta1)); ...
                    -sin(theta1(itheta1)) cos(theta1(itheta1));   ];
            W([i j],:)= c * W([i j],:);
            param.fixedtree=0;
            [score,details]=contrast_tca(xnorm,W,param);
            if (itheta1<=2) | (itheta1>=searchgrid),
                finished(i,j)=1;
                if sum(finished(:))==m*(m-1)/2, return; end
            else
                finished=zeros(m);
            end
        end
    end
 if detail,   fprintf('sweep %d, oldscore=%f newscore=%f\n',k,reallyoldscore,score); end
    
end