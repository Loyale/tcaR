function [val,par] = sphere_tca(theta,type,par,xnorm,W,param,details,i,j,invD,theta10,theta20)

par=[];
theta1=theta(1);
theta2=theta(2);


c = [cos(theta1+theta10) sin(theta1+theta10); ...
        cos(theta2+theta20) sin(theta2+theta20)];
a = c * invD;
W([i j],:)=a * W([i j],:);
val=contrast_tca_multupdate(xnorm,W,param,[i j],details);         