function [cost, Tf] = segment_cost(from_node, to_point, dim)

if dim==2
        
    Tf = roots(DI_time(from_node(1), from_node(2), from_node(3), from_node(4), to_point(1), to_point(2), to_point(3), to_point(4)));
    Tf = Tf(abs(imag(Tf)) < 0.0001);
    Tf = min(Tf(Tf >= 0));

    cost = DI_cost(Tf, from_node(1), from_node(2), from_node(3), from_node(4), to_point(1), to_point(2), to_point(3), to_point(4));
    
elseif dim == 3
    
    Tf = roots(DI3d_time(from_node(1), from_node(2), from_node(3), from_node(4), from_node(5), from_node(6), to_point(1), to_point(2), to_point(3), to_point(4), to_point(5), to_point(6)));
        Tf = Tf(abs(imag(Tf)) < 0.0001);
        Tf = min(Tf(Tf >= 0));

    cost = DI3d_cost(Tf, from_node(1), from_node(2), from_node(3), from_node(4), from_node(5), from_node(6), to_point(1), to_point(2), to_point(3), to_point(4), to_point(5), to_point(6));
    
end
    
end