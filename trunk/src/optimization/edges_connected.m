function list_edges = edges_connected(triang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Given a triangulation, calculate the list of the vertices that are connected
%% to each of its vertices (triang.X), i.e, it gives all the edges.
%%
%% Inputs:  1. triang -> triangulation class X & Triangulation
%% Outputs: 1. list_edges -> cell of size equals to the number of points in 
%%                           the given triang
%%
%% NOTE: some edges are repeted, because here the order matters
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

triang_edges = edges(triang);
list_edges = cell(size(triang.X,1),1);

for i = 1:size(triang.X,1)
    
    [r_tmp,c_tmp] = find(triang_edges == i);
    list_edges{i} = [];
    
    for j = 1:length(r_tmp)
        list_edges{i} = [list_edges{i} triang_edges(r_tmp(j),c_tmp(j) + (-1)^(c_tmp(j)-1))];%diag(triang_edges(r_tmp,c_tmp));
    end
    
end

