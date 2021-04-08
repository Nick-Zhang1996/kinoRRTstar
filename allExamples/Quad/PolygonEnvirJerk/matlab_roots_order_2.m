% solve 2nd order polynomial, give smallest positive solution
% takes coefficients, high to low order
function [res] = matlab_roots_order_2(p1,p2,p3) %codegen
    poly = [p1,p2,p3];
    Tf = roots(poly);
    mask = abs(imag(Tf)) < 0.0001;
    Tf = Tf(mask);
    Tf = Tf(Tf >= 0);
    res = Tf(1);
end