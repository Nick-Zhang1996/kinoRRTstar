function cost = DI3d_costFreeVel(t_s,x01,x02,x03,x04,x05,x06,x11,x12,x13)
    %sprintf("DI3d_costFreeVel(%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f)",t_s,x01,x02,x03,x04,x05,x06,x11,x12,x13)
    cost = t_s + (3*x04^2 + 3*x05^2 + 3*x06^2)/t_s + (6*x01*x04 + 6*x02*x05 + 6*x03*x06 - 6*x04*x11 ...
        - 6*x05*x12 - 6*x06*x13)/t_s^2 + (3*x01^2 - 6*x01*x11 + 3*x02^2 - 6*x02*x12 + 3*x03^2 - 6*x03*x13 + 3*x11^2 + 3*x12^2 + 3*x13^2)/t_s^3;
end