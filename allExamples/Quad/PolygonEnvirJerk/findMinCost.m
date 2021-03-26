function [cost] = findMinCost(tree,end_node,dim,state_dim)

    idx = tree(:,state_dim+1)==1;
    connectingNodes = tree(idx,:);
    
    temp = size(connectingNodes,1);
    PathCost = zeros(temp,1);
    
    for k = 1:temp
        [lcost, ~] = segment_cost(connectingNodes(k, :), end_node, dim);
        PathCost(k) = connectingNodes(k, state_dim+2) + real(lcost);
    end           
    
    % find minimum cost last node
    cost = min(PathCost);
    cost = real(cost);
end