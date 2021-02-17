%calculate the cost from root to to_point through parent from_node
function [cost, Tf] = cost_np(from_node, to_point, dim, Tf)

if dim == 2

%     if ~exist('Tf','var') || isempty(Tf)                
%         Tf = roots(obj.eval_arrival_internal(from_node(1), from_node(2), from_node(3), from_node(4), to_point(1), to_point(2), to_point(3), to_point(4)));
%         Tf = Tf(imag(Tf) == 0);
%         Tf = min(Tf(Tf >= 0));
%     end

    cost = cost_eval(Tf, from_node(1), from_node(2), from_node(3), from_node(4), to_point(1), to_point(2), to_point(3), to_point(4));
    cost = from_node(:, 2 * dim + 2) + cost;

    
elseif dim == 3
    cost = 100;
    Tf = 0;    
end

end