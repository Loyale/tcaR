function [components,var_comp_matrix]=connected_components(G);
EG=expm(G);
N=length(G);
not_in_a_component=ones(1,N);

k=1;
while (k<=N)
    ind=find(not_in_a_component);
    if isempty(ind), break; end
    var_comp_matrix(:,k)=(EG(:,ind(1))>0);
    components{k}=find(var_comp_matrix(:,k));
    not_in_a_component(components{k})=0;
    k=k+1;
end
