function [mypath, cost] = findMinimumPath(tree,end_node,dim)
    %NOTE
    coder.varsize("mypath");
    % find nodes that connect to end_node
%     connectingNodes = [];
%     for i=1:size(tree,1)
%         if tree(i,2*dim+1)==1
%             tree(i,2*dim+2)=tree(i,2*dim+2)+segment_cost(tree(i,:),end_node(1:2*dim),dim);
%             connectingNodes = [connectingNodes ; tree(i,:)];
%         end
%     end
    idx = tree(:,2*dim+1)==1;
    connectingNodes = tree(idx,:);
    %NOTE
    temp = size(connectingNodes);
    for k = 1:temp(1)
        if  norm(connectingNodes(k, 1 : dim) - end_node(1 : dim)) >= 0.2
            [lcost, ~] = segment_cost(connectingNodes(k, :), end_node, dim);
            connectingNodes(k, 2*dim+2) = connectingNodes(k, 2*dim+2) + lcost;
        end
    end           
    
    % find minimum cost last node
    [cost,idx] = min(connectingNodes(:,2*dim+2));
    
    % construct lowest cost path
    if  norm(connectingNodes(idx, 1 : dim) - end_node(1 : dim)) >= 0.2
        [lcost, ~] = segment_cost(connectingNodes(idx, :), end_node, dim);
        connectingNodes(idx, 2*dim+2) = connectingNodes(idx, 2*dim+2) - lcost;
        mypath = [connectingNodes(idx,:); end_node];
        
        parent_node = connectingNodes(idx,2*dim+3);
        while parent_node>=1
            mypath = [tree(parent_node,:); mypath];
            parent_node = tree(parent_node,2*dim+3);       
        end
    else
        mypath = [connectingNodes(idx,:)];
        parent_node = connectingNodes(idx,2*dim+3);
        while parent_node>=1
            mypath = [tree(parent_node,:); mypath];
            parent_node = tree(parent_node,2*dim+3);       
        end    
    end

end