function [score,details]=contrast_tca_kgv(x,W,tparam);


[m N]=size(x);
if (nargout>1), detail=1; else detail=0; end
x=W*x;
sigma=x*x'/N; % assumes zero-mean !!!!!!!!!!!11


lambdaC=tparam.lambdaC;
lambdaG=tparam.lambdaG;
lambdaT=tparam.lambdaT;

kkappa=tparam.kappa;
keta=tparam.eta;

Rkappa=[];
sizes=[];

for i=1:m
    % cholesky decomposition using a MEX-file
    switch (tparam.kernel)
        case 'hermite'
            [G,Pvec] =chol_hermite(x(i,:),tparam.sigma,tparam.p,N*tparam.eta); 
        case 'gaussian'
            [G,Pvec] =chol_gauss(x(i,:)/tparam.sigma,1,N*tparam.eta); 
        case 'lineargaussian'
            [G,Pvec] =chol_gauss(x(i,:),tparam.sigma,tparam.lambdaKG,N*tparam.eta); 
        case 'poly'
            [G,Pvec] =chol_poly(x(i,:),tparam.r,tparam.s,tparam.d,N*tparam.eta); 
    end
    
    [a,Pvec]=sort(Pvec);
    G=centerpartial(G(Pvec,:));
    
    % regularization (see paper for details)
    [A,D]=eig(G'*G);
    D=diag(D);
    indexes=find(D>=N*keta & isreal(D)); %removes small eigenvalues
    [newinds,order]=sort(D(indexes));
    order=flipud(order);
    neig=length(indexes);
    indexes=indexes(order(1:neig));  
    if (isempty(indexes)), indexes=[1]; end
    D=D(indexes);
    V=G*(A(:,indexes)*diag(sqrt(1./(D))));
    Us{i}=V;
    Lambdas{i}=D;
    Dr=D;
    for j=1:length(D)
        Dr(j)=D(j)/(N*kkappa+D(j));
    end
    Drs{i}=Dr;
    sizes(i)=size(Drs{i},1);
end




% calculate Rkappa
Rkappa=eye(sum(sizes));
starts=cumsum([1 sizes]);
starts(m+1)=[];
for i=2:m
    for j=1:(i-1)
        newbottom=diag(Drs{i})*(Us{i}'*Us{j})*diag(Drs{j});
        Rkappa(starts(i):starts(i)+sizes(i)-1,starts(j):starts(j)+sizes(j)-1)=newbottom;
        Rkappa(starts(j):starts(j)+sizes(j)-1,starts(i):starts(i)+sizes(i)-1)=newbottom';
    end
end
if (detail)
    % outputs details
    details.Us=Us;
    details.Lambdas=Lambdas;
    details.Drs=Drs;
    details.Rkappa=Rkappa;
    details.sizes=sizes;
    details.starts=starts;
end

% MUTUAL INFORMATION
J=-.5*log(det(Rkappa));
details.mutinf=J;



weights=zeros(m);
% split, whether we need to calculate the tree or not
% compute all pairwise mutual informations
weights=zeros(m);
for i=1:m-1
    for j=i+1:m
        Rkappap1=Rkappa(starts(i):starts(i)+sizes(i)-1,starts(j):starts(j)+sizes(j)-1);
        Rkappasm=[eye(sizes(i)) Rkappap1; Rkappap1' eye(sizes(j))];
        weights(i,j)=-.5*log(det(Rkappasm));
        corr=sigma(i,j)/sigma(i,i)^.5/sigma(j,j)^.5;
        
        weights(i,j)=weights(i,j)-lambdaG*.5*log(1-corr^2)+ ...
            +lambdaC*.5*log(1+1e-14-corr^2);
        weights(j,i)=weights(i,j);
    end
end
details.weights=weights;
if isfield(tparam,'edgeprior'), edgeprior=tparam.edgeprior; else edgeprior=0; end
[tree,weightsum] = minimum_nonspanning_tree_edges(edgeprior-weights,[],tparam.maxedges);
    J=J+lambdaT*weightsum;

details.tree=tree;
score=J;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G2=centerpartial(G1)
% CENTERPARTIAL - Center a gram matrix of the form K=G*G'

[N,NG]=size(G1);
G2 = G1 - repmat(mean(G1,1),N,1);







