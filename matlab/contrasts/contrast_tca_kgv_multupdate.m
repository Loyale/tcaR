function      [score,details]=contrast_tca_kgv_multupdate(x,W,param,js,details);

if (nargout>1), detail=1; else detail=0; end

[m N]=size(x);
x=W*x;
sigma=x*x'/N;

kkappa=param.kappa;
keta=param.eta;
lambdaG=param.lambdaG;
lambdaC=param.lambdaC;
lambdaT=param.lambdaT;

Rkappa=details.Rkappa;
Us=details.Us;
Lambdas=details.Lambdas;
Drs=details.Drs;
sizes=details.sizes;
oldstarts=details.starts;
oldsizes=sizes;
weights=details.weights;

%redo the 2 cholesky decompositions
for k=js
    switch (param.kernel)
        case 'hermite'
            [G,Pvec] =chol_hermite(x(k,:),param.sigma,param.p,N*param.eta); 
        case 'gaussian'
            [G,Pvec] =chol_gauss(x(k,:)/param.sigma,1,N*param.eta); 
        case 'lineargaussian'
            [G,Pvec] =chol_gauss(x(k,:),param.sigma,param.lambdaKG,N*param.eta); 
            
        case 'poly'
            [G,Pvec] =chol_poly(x(k,:),param.r,param.s,param.d,N*param.eta); 
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
    Us{k}=V;
    Lambdas{k}=D;
    Dr=D;
    for j2=1:length(D)
        Dr(j2)=D(j2)/(N*kkappa+D(j2));
    end
    Drs{k}=Dr;
    sizes(k)=size(Drs{k},1);
end
starts=cumsum([1 sizes]);
starts(m+1)=[];

newRkappa=eye(sum(sizes));

for i=2:m
    for j=1:i-1
        if ( ismember(i,js) | ismember(j,js) )
            newbottom=diag(Drs{i})*(Us{i}'*Us{j})*diag(Drs{j});
            newRkappa(starts(i):starts(i)+sizes(i)-1,starts(j):starts(j)+sizes(j)-1)=newbottom;
            newRkappa(starts(j):starts(j)+sizes(j)-1,starts(i):starts(i)+sizes(i)-1)=newbottom';
        else
            newbottom= Rkappa(oldstarts(i):oldstarts(i)+oldsizes(i)-1,oldstarts(j):oldstarts(j)+oldsizes(j)-1);
            newRkappa(starts(i):starts(i)+sizes(i)-1,starts(j):starts(j)+sizes(j)-1)=newbottom;
            newRkappa(starts(j):starts(j)+sizes(j)-1,starts(i):starts(i)+sizes(i)-1)=newbottom';
        end
    end
end

Rkappa=newRkappa;
clear newRkappa;
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

if param.fixedtree
    tree=param.tree;
    weightsum=0;
    for k=1:size(param.tree,2);
        weightsum = weightsum - weights(param.tree(1,k),param.tree(2,k));
    end
else
if isfield(param,'edgeprior'), edgeprior=param.edgeprior; else edgeprior=0; end
[tree,weightsum] = minimum_nonspanning_tree_edges(edgeprior-weights,[],param.maxedges);
end

J=J+lambdaT*weightsum;

details.tree=tree;
score=J;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G2=centerpartial(G1)
% CENTERPARTIAL - Center a gram matrix of the form K=G*G'

[N,NG]=size(G1);
G2 = G1 - repmat(mean(G1,1),N,1);


