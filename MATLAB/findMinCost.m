function [min_cost] = findMinCost(tree,end_node,dim)
    % find nodes that connect to end_node
%     connectingNodes = [];
%     for i=1:size(tree,1)
%         if tree(i,2*dim+1)==1
%             tree(i,2*dim+2)=tree(i,2*dim+2)+segment_cost(tree(i,:),end_node(1:2*dim),dim);
%             connectingNodes = [connectingNodes ; tree(i,:)];
%         end
%     end
    close_to_goal_idx = tree(:,2*dim+1)==1;
    %fprintf("close to goal index %.1f \n",sum(nonzeros(close_to_goal_idx)));
    connectingNodes = tree(close_to_goal_idx,:);
    %DEBUG
    i=1;
    %fprintf("connecting Nodes  : %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, \n",connectingNodes(i,1),connectingNodes(i,2),connectingNodes(i,3),connectingNodes(i,4),connectingNodes(i,5),connectingNodes(i,6),connectingNodes(i,7),connectingNodes(i,8));
    min_cost = 10000;
    %fprintf("connecting Nodes size %.0f",size(connectingNodes,1));
    for k = 1:size(connectingNodes,1)
        if  norm(connectingNodes(k, 1 : dim) - end_node(1 : dim)) >= 0.2

            [lcost, ~] = segment_cost(connectingNodes(k, :), end_node, dim);
            %fprintf("evaluate node cost, local to end = %.2f \n",lcost);
            total_cost = connectingNodes(k, 2*dim+2) + real(lcost);
            %fprintf("evaluate node cost, total = %.2f \n",total_cost);
            if total_cost < min_cost
                min_cost = total_cost;
            end
        end
    end           
    

end