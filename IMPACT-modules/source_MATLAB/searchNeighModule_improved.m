function Neigh_all = searchNeighModule_improved(A,nodes)
% function Neigh_all = searchNeighModule_improved ( A , nodes )
% A :         Adjacency matrix (note -- it hase to be the complete adj matrix and not only the upper or lower triangular part
% nodes:      current set of query nodes for finding their neigh
%

if length(nodes)==1
    ind_module=1;
elseif length(nodes)>1
    ind_module=2;
end
    

    neigh_all=[];

    for i=ind_module:length(nodes)

        curr_node=nodes(i);
        curr_neigh1=find(A(curr_node,:)==1);
        curr_neigh2=find(A(:,curr_node)==1);
        curr_neigh=unique(union(curr_neigh1,curr_neigh2));
        if size(curr_neigh,1)>size(curr_neigh,2)
            curr_neigh=curr_neigh';
        end
        
        neigh_all=[neigh_all curr_neigh];
    end

Neigh_all=unique(neigh_all);