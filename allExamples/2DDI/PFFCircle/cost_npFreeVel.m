%calculate the cost from root to to_point through parent from_node
function [cost, Tf] = cost_npFreeVel(from_node, to_point, dim, Tf)

if dim == 2

    cost = DI_costFreeVel(Tf, from_node(1), from_node(2), from_node(3), from_node(4), to_point(1), to_point(2));
    %fprintf("DI_costFreeVel %.1f\n", cost);
%     fprintf("    from node ");
%     for i = 1:6
%         fprintf(" %.2f, ",from_node(i));
%     end
%     fprintf("\n");
    
    %fprintf("    parent cost %.1f\n", from_node( 2 * dim + 2));
    % NOTE do not use broadcasting
    cost = from_node( 2 * dim + 2) + cost;
    
elseif dim == 3
    cost = DI3d_costFreeVel(Tf, from_node(1), from_node(2), from_node(3), from_node(4), from_node(5), from_node(6), to_point(1), to_point(2), to_point(3));
    cost = from_node( 2 * dim + 2) + cost;
end

end