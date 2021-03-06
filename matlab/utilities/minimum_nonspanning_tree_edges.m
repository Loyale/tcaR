function [A,weightsum] = minimum_nonspanning_tree_edges(C1, C2, maxedges)
%
% Find the minimum spanning tree using Prim's algorithm.
% C1(i,j) is the primary cost of connecting i to j.
% C2(i,j) is the (optional) secondary cost of connecting i to j, used to break ties.
% We assume that absent edges have 0 cost.
% To find the maximum spanning tree, used -1*C.
% See Aho, Hopcroft & Ullman 1983, "Data structures and algorithms", p 237.

% Prim's is O(V^2). Kruskal's algorithm is O(E log E) and hence is more efficient
% for sparse graphs, but is implemented in terms of a priority queue.

% We partition the nodes into those in U and those not in U.
% closest(i) is the vertex in U that is closest to i in V-U.
% lowcost(i) is the cost of the edge (i, closest(i)), or infinity is i has been used.
% In Aho, they say C(i,j) should be "some appropriate large value" if the edge is missing.
% We set it to infinity.
% However, since lowcost is initialized from C, we must distinguish absent edges from used nodes.
weights=[];
weightsum = 0;
n = length(C1);
if nargin==1, C2 = zeros(n);  maxedges=n-1; end
if nargin==2, maxedges=n-1; end
if isempty(C2), C2 = zeros(n); end
% A = zeros(n);
A = [];

closest = ones(1,n);
used = zeros(1,n); % contains the members of U
used(1) = 1; % start with node 1
% C1(find(C1==0))=inf;
% C2(find(C2==0))=inf;
C1=setdiag(C1,inf*ones(n,1));
C2=setdiag(C2,inf*ones(n,1));

lowcost1 = C1(1,:);
lowcost2 = C2(1,:);

for i=2:n
  ks = find(lowcost1==min(lowcost1));
  k = ks(argmin(lowcost2(ks)));
%   A(k, closest(k)) = 1;
%   A(closest(k), k) = 1;
A = [ A [k; closest(k)]];
weightsum = weightsum + C1(k,closest(k));
weights = [weights C1(k,closest(k))];
  lowcost1(k) = inf;
  lowcost2(k) = inf;
  used(k) = 1;
  NU = find(used==0);
  for ji=1:length(NU)
    for j=NU(ji)
      if C1(k,j) < lowcost1(j)
	lowcost1(j) = C1(k,j);
	lowcost2(j) = C2(k,j);
	closest(j) = k;
      end
    end
  end
end
ind=find(weights<0);
A=A(:,ind);
weights=weights(ind);
weightsum=sum(weights);    
[a,b]=sort(weights);
if length(b)>=maxedges,
A=A(:,b(1:maxedges));
weightsum=sum(weights(b(1:maxedges)));    
end

