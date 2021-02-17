function cost =  cost_eval(t_s,x01,x02,x03,x04,x11,x12,x13,x14)
    cost = t_s + (4*x03^2)/t_s + (12*x01^2)/t_s^3 + (4*x04^2)/t_s + (12*x02^2)/t_s^3 +...
    (4*x13^2)/t_s + (12*x11^2)/t_s^3 + (4*x14^2)/t_s + (12*x12^2)/t_s^3 +...
    (12*x01*x03)/t_s^2 + (12*x02*x04)/t_s^2 - (24*x01*x11)/t_s^3 + (12*x01*x13)/t_s^2 -...
    (12*x03*x11)/t_s^2 + (4*x03*x13)/t_s - (24*x02*x12)/t_s^3 + (12*x02*x14)/t_s^2 -...
    (12*x04*x12)/t_s^2 + (4*x04*x14)/t_s - (12*x11*x13)/t_s^2 - (12*x12*x14)/t_s^2;
    %DEBUG
    %fprintf("cost_eval %.2f, %.2f, %.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f, \n",t_s,x01,x02,x03,x04,x11,x12,x13,x14);
end