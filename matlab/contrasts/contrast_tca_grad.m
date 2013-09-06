function [J0,gradJs]=contrast_tca_grad(x,W,param)

m=size(W,1);
gradJs=cell(1,m);
[J0,details]=contrast_tca(x,W,param);
s=W'*x;
fixedind=param.fixedind;
variableind=setdiff(1:m,param.fixedind);


for i=variableind
    w0=W(:,i);
    wc=null(w0');
    WTgradF=zeros(m,1);
    dr=0.001;
    for j=2:m
        wdr=w0*cos(dr)+wc(:,j-1)*sin(dr);
        Wnew=W;
        Wnew(:,i)=wdr;
        J=contrast_tca_update(x,Wnew,param,i,details);
        WTgradF=WTgradF+(J-J0)/dr*wc(:,j-1);
    end
   % WTgradF
    gradJs{i}=WTgradF;
    if (1)
        WTgradF=zeros(m,1);
        dr=-0.001;
        for j=2:m
            wdr=w0*cos(dr)+wc(:,j-1)*sin(dr);
            Wnew=W;
            Wnew(:,i)=wdr;
            J=contrast_tca_update(x,Wnew,param,i,details);
            WTgradF=WTgradF+(J-J0)/dr*wc(:,j-1);
        end
      %  WTgradF
      %  pause
        gradJs{i}=(WTgradF+gradJs{i})/2;
    end
end   