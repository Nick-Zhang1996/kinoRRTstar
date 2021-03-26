function cost = quadf_cost(t_s,x01, x02, x03, x04, x05, x06, x07, x08, x09, x11, x12, x13, x14, x15, x16, x17, x18, x19)
cost =t_s + (18*x01^2)/(125*t_s^5) + (24*x04^2)/(625*t_s^3) + (18*x02^2)/(125*t_s^5) + (9*x07^2)/(5000*t_s) + (24*x05^2)/(625*t_s^3) + (18*x03^2)/(125*t_s^5) + (9*x08^2)/(5000*t_s) + (24*x06^2)/(625*t_s^3) + (9*x09^2)/(5000*t_s) + (18*x11^2)/(125*t_s^5) + (24*x14^2)/(625*t_s^3) + (18*x12^2)/(125*t_s^5) + (9*x17^2)/(5000*t_s) + (24*x15^2)/(625*t_s^3) + (18*x13^2)/(125*t_s^5) + (9*x18^2)/(5000*t_s) + (24*x16^2)/(625*t_s^3) + (9*x19^2)/(5000*t_s) + (18*x01*x04)/(125*t_s^4) + (3*x01*x07)/(125*t_s^3) + (18*x02*x05)/(125*t_s^4) + (9*x04*x07)/(625*t_s^2) + (3*x02*x08)/(125*t_s^3) + (18*x03*x06)/(125*t_s^4) + (9*x05*x08)/(625*t_s^2) + (3*x03*x09)/(125*t_s^3) + (9*x06*x09)/(625*t_s^2) - (36*x01*x11)/(125*t_s^5) + (18*x01*x14)/(125*t_s^4) - (18*x04*x11)/(125*t_s^4) - (36*x02*x12)/(125*t_s^5) - (3*x01*x17)/(125*t_s^3) + (42*x04*x14)/(625*t_s^3) - (3*x07*x11)/(125*t_s^3) + (18*x02*x15)/(125*t_s^4) - (18*x05*x12)/(125*t_s^4) - (36*x03*x13)/(125*t_s^5) - (6*x04*x17)/(625*t_s^2) + (6*x07*x14)/(625*t_s^2) - (3*x02*x18)/(125*t_s^3) + (42*x05*x15)/(625*t_s^3) - (3*x08*x12)/(125*t_s^3) + (18*x03*x16)/(125*t_s^4) - (18*x06*x13)/(125*t_s^4) - (3*x07*x17)/(2500*t_s) - (6*x05*x18)/(625*t_s^2) + (6*x08*x15)/(625*t_s^2) - (3*x03*x19)/(125*t_s^3) + (42*x06*x16)/(625*t_s^3) - (3*x09*x13)/(125*t_s^3) - (3*x08*x18)/(2500*t_s) - (6*x06*x19)/(625*t_s^2) + (6*x09*x16)/(625*t_s^2) - (3*x09*x19)/(2500*t_s) - (18*x11*x14)/(125*t_s^4) + (3*x11*x17)/(125*t_s^3) - (18*x12*x15)/(125*t_s^4) - (9*x14*x17)/(625*t_s^2) + (3*x12*x18)/(125*t_s^3) - (18*x13*x16)/(125*t_s^4) - (9*x15*x18)/(625*t_s^2) + (3*x13*x19)/(125*t_s^3) - (9*x16*x19)/(625*t_s^2);
cost = real(cost);
end