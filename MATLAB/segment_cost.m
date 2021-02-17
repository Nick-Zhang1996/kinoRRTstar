function [cost, Tf] = segment_cost(from_node, to_point, dim)

if dim==2
        
        Tf = roots(opt_time(from_node(1), from_node(2), from_node(3), from_node(4), to_point(1), to_point(2), to_point(3), to_point(4)));
        Tf = Tf(abs(imag(Tf)) < 0.0001);
        %NOTE
        Tf = min(Tf(Tf >= 0));

    cost = cost_eval(Tf, from_node(1), from_node(2), from_node(3), from_node(4), to_point(1), to_point(2), to_point(3), to_point(4));

elseif dim == 3
    cost = 100;
    Tf = 0;   
end
    
end