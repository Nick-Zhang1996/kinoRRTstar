function [path, cost] = findMinimumPath(tree,end_node,dim,state_dim)

    coder.varsize("path");

    idx = tree(:,state_dim+1)==1;
    connectingNodes = tree(idx,:);
    
    temp = size(connectingNodes,1);
    PathCost = zeros(temp,1);
    
    for k = 1:temp
        [lcost, ~] = segment_cost(connectingNodes(k, :), end_node, dim);
        PathCost(k) = connectingNodes(k, state_dim+2) + lcost;
    end           
    
    % find minimum cost last node
    [cost,idx] = min(PathCost);
    
    % construct lowest cost path        
    path = [connectingNodes(idx,:); end_node];
    parent_node = connectingNodes(idx,state_dim+3);
    while parent_node>=1
        path = [tree(parent_node,:); path];
        parent_node = tree(parent_node,state_dim+3);       
    end


end