function score=orth_tca(theta,xnorm, W,param,details,i,j);

c = [cos(theta) sin(theta); ...
                    -sin(theta) cos(theta);   ];
            W([i j],:)= c * W([i j],:);
  score=contrast_tca_multupdate(xnorm,W,param,[i j],details);