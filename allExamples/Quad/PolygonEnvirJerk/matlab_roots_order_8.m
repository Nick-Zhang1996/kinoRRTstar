% solve 8th order polynomial, give smallest positive solution
% takes coefficients, high to low order
function [res] = matlab_roots_order_8(p1,p2,p3,p4,p5,p6,p7,p8,p9)
    poly = [p1,p2,p3,p4,p5,p6,p7,p8,p9];
    Tf = roots(poly);
    Tf = Tf(abs(imag(Tf)) < 0.0001);
    Tf = Tf(Tf >= 0);
    res = Tf(1);
end