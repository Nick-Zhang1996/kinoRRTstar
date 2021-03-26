function [cost, Tf] = segment_cost(from_node, to_point, dim)

if dim==2
        
%     Tf = norm(from_node(1 : 4) - to_point(1 : 4))/0.8;
    Tf = roots(DI_time(from_node(1), from_node(2), from_node(3), from_node(4), to_point(1), to_point(2), to_point(3), to_point(4)));
    Tf = Tf(abs(imag(Tf)) < 0.0001);
    Tf = min(Tf(Tf >= 0));

    cost = DI_cost(Tf, from_node(1), from_node(2), from_node(3), from_node(4), to_point(1), to_point(2), to_point(3), to_point(4));
    
elseif dim == 3
    
    Tf = roots(quadf_time(from_node(1), from_node(2), from_node(3), from_node(4), from_node(5), from_node(6), from_node(7), from_node(8), from_node(9), to_point(1), to_point(2), to_point(3), to_point(4), to_point(5), to_point(6), to_point(7), to_point(8), to_point(9)));
    Tf = Tf(abs(imag(Tf)) < 0.0001);
    Tf = Tf(Tf >= 0);
    tf = Tf(1);
    cost = quadf_cost(Tf(1), from_node(1), from_node(2), from_node(3), from_node(4), from_node(5), from_node(6), from_node(7), from_node(8), from_node(9), to_point(1), to_point(2), to_point(3), to_point(4), to_point(5), to_point(6), to_point(7), to_point(8), to_point(9));
    for i = 2 : length(Tf)
        costnow = quadf_cost(Tf(i), from_node(1), from_node(2), from_node(3), from_node(4), from_node(5), from_node(6), from_node(7), from_node(8), from_node(9), to_point(1), to_point(2), to_point(3), to_point(4), to_point(5), to_point(6), to_point(7), to_point(8), to_point(9));
        if costnow<cost
            tf = Tf(i);
            cost = costnow;
        end
    end
    Tf = tf;
    
%     Tf = norm(from_node(1 : 3) - to_point(1 : 3))/4;
    
    cost = quadf_cost(Tf(1), from_node(1), from_node(2), from_node(3), from_node(4), from_node(5), from_node(6), from_node(7), from_node(8), from_node(9), to_point(1), to_point(2), to_point(3), to_point(4), to_point(5), to_point(6), to_point(7), to_point(8), to_point(9));
    
end
    
end