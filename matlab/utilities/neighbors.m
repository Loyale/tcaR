function ns = neighbors(adj_mat, i)

% Finds the neighbors of a node in a graph

  ns = find(adj_mat(i,:));

