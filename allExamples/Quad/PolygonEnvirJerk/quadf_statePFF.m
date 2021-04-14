function states = quadf_statePFF(t,t_s,x01,x02,x03,x04,x05,x06,x07,x08,x09,x11,x12,x13)
states = [x01 + t*x04 + (t^2*x07)/2 + (625*t^3*x01)/(500*t_s^3 + 9*t_s^2) - (15*t^3*x01)/(500*t_s^4 + 9*t_s^3) + (15*t^4*x01)/(2*(500*t_s^5 + 9*t_s^4)) - (3*t^5*x01)/(2*(500*t_s^6 + 9*t_s^5)) - (625*t^3*x11)/(500*t_s^3 + 9*t_s^2) + (15*t^3*x11)/(500*t_s^4 + 9*t_s^3) - (15*t^4*x11)/(2*(500*t_s^5 + 9*t_s^4)) + (3*t^5*x11)/(2*(500*t_s^6 + 9*t_s^5)) - (1875*t^3*t_s*x01)/(500*t_s^4 + 9*t_s^3) + (625*t^3*t_s*x04)/(500*t_s^3 + 9*t_s^2) - (15*t^3*t_s*x04)/(500*t_s^4 + 9*t_s^3) + (1875*t^4*t_s*x01)/(2*(500*t_s^5 + 9*t_s^4)) + (15*t^4*t_s*x04)/(2*(500*t_s^5 + 9*t_s^4)) - (375*t^5*t_s*x01)/(2*(500*t_s^6 + 9*t_s^5)) - (3*t^5*t_s*x04)/(2*(500*t_s^6 + 9*t_s^5)) + (1875*t^3*t_s*x11)/(500*t_s^4 + 9*t_s^3) - (1875*t^4*t_s*x11)/(2*(500*t_s^5 + 9*t_s^4)) + (375*t^5*t_s*x11)/(2*(500*t_s^6 + 9*t_s^5)) - (1875*t^3*t_s^2*x04)/(500*t_s^4 + 9*t_s^3) + (125*t^3*t_s^2*x07)/(500*t_s^3 + 9*t_s^2) - (15*t^3*t_s^2*x07)/(2*(500*t_s^4 + 9*t_s^3)) + (1875*t^4*t_s^2*x04)/(2*(500*t_s^5 + 9*t_s^4)) - (625*t^3*t_s^3*x07)/(500*t_s^4 + 9*t_s^3) + (15*t^4*t_s^2*x07)/(4*(500*t_s^5 + 9*t_s^4)) - (375*t^5*t_s^2*x04)/(2*(500*t_s^6 + 9*t_s^5)) + (625*t^4*t_s^3*x07)/(2*(500*t_s^5 + 9*t_s^4)) - (3*t^5*t_s^2*x07)/(4*(500*t_s^6 + 9*t_s^5)) - (125*t^5*t_s^3*x07)/(2*(500*t_s^6 + 9*t_s^5));
 x02 + t*x05 + (t^2*x08)/2 + (625*t^3*x02)/(500*t_s^3 + 9*t_s^2) - (15*t^3*x02)/(500*t_s^4 + 9*t_s^3) + (15*t^4*x02)/(2*(500*t_s^5 + 9*t_s^4)) - (3*t^5*x02)/(2*(500*t_s^6 + 9*t_s^5)) - (625*t^3*x12)/(500*t_s^3 + 9*t_s^2) + (15*t^3*x12)/(500*t_s^4 + 9*t_s^3) - (15*t^4*x12)/(2*(500*t_s^5 + 9*t_s^4)) + (3*t^5*x12)/(2*(500*t_s^6 + 9*t_s^5)) - (1875*t^3*t_s*x02)/(500*t_s^4 + 9*t_s^3) + (625*t^3*t_s*x05)/(500*t_s^3 + 9*t_s^2) - (15*t^3*t_s*x05)/(500*t_s^4 + 9*t_s^3) + (1875*t^4*t_s*x02)/(2*(500*t_s^5 + 9*t_s^4)) + (15*t^4*t_s*x05)/(2*(500*t_s^5 + 9*t_s^4)) - (375*t^5*t_s*x02)/(2*(500*t_s^6 + 9*t_s^5)) - (3*t^5*t_s*x05)/(2*(500*t_s^6 + 9*t_s^5)) + (1875*t^3*t_s*x12)/(500*t_s^4 + 9*t_s^3) - (1875*t^4*t_s*x12)/(2*(500*t_s^5 + 9*t_s^4)) + (375*t^5*t_s*x12)/(2*(500*t_s^6 + 9*t_s^5)) - (1875*t^3*t_s^2*x05)/(500*t_s^4 + 9*t_s^3) + (125*t^3*t_s^2*x08)/(500*t_s^3 + 9*t_s^2) - (15*t^3*t_s^2*x08)/(2*(500*t_s^4 + 9*t_s^3)) + (1875*t^4*t_s^2*x05)/(2*(500*t_s^5 + 9*t_s^4)) - (625*t^3*t_s^3*x08)/(500*t_s^4 + 9*t_s^3) + (15*t^4*t_s^2*x08)/(4*(500*t_s^5 + 9*t_s^4)) - (375*t^5*t_s^2*x05)/(2*(500*t_s^6 + 9*t_s^5)) + (625*t^4*t_s^3*x08)/(2*(500*t_s^5 + 9*t_s^4)) - (3*t^5*t_s^2*x08)/(4*(500*t_s^6 + 9*t_s^5)) - (125*t^5*t_s^3*x08)/(2*(500*t_s^6 + 9*t_s^5));
 x03 + t*x06 + (t^2*x09)/2 + (625*t^3*x03)/(500*t_s^3 + 9*t_s^2) - (15*t^3*x03)/(500*t_s^4 + 9*t_s^3) + (15*t^4*x03)/(2*(500*t_s^5 + 9*t_s^4)) - (3*t^5*x03)/(2*(500*t_s^6 + 9*t_s^5)) - (625*t^3*x13)/(500*t_s^3 + 9*t_s^2) + (15*t^3*x13)/(500*t_s^4 + 9*t_s^3) - (15*t^4*x13)/(2*(500*t_s^5 + 9*t_s^4)) + (3*t^5*x13)/(2*(500*t_s^6 + 9*t_s^5)) - (1875*t^3*t_s*x03)/(500*t_s^4 + 9*t_s^3) + (625*t^3*t_s*x06)/(500*t_s^3 + 9*t_s^2) - (15*t^3*t_s*x06)/(500*t_s^4 + 9*t_s^3) + (1875*t^4*t_s*x03)/(2*(500*t_s^5 + 9*t_s^4)) + (15*t^4*t_s*x06)/(2*(500*t_s^5 + 9*t_s^4)) - (375*t^5*t_s*x03)/(2*(500*t_s^6 + 9*t_s^5)) - (3*t^5*t_s*x06)/(2*(500*t_s^6 + 9*t_s^5)) + (1875*t^3*t_s*x13)/(500*t_s^4 + 9*t_s^3) - (1875*t^4*t_s*x13)/(2*(500*t_s^5 + 9*t_s^4)) + (375*t^5*t_s*x13)/(2*(500*t_s^6 + 9*t_s^5)) - (1875*t^3*t_s^2*x06)/(500*t_s^4 + 9*t_s^3) + (125*t^3*t_s^2*x09)/(500*t_s^3 + 9*t_s^2) - (15*t^3*t_s^2*x09)/(2*(500*t_s^4 + 9*t_s^3)) + (1875*t^4*t_s^2*x06)/(2*(500*t_s^5 + 9*t_s^4)) - (625*t^3*t_s^3*x09)/(500*t_s^4 + 9*t_s^3) + (15*t^4*t_s^2*x09)/(4*(500*t_s^5 + 9*t_s^4)) - (375*t^5*t_s^2*x06)/(2*(500*t_s^6 + 9*t_s^5)) + (625*t^4*t_s^3*x09)/(2*(500*t_s^5 + 9*t_s^4)) - (3*t^5*t_s^2*x09)/(4*(500*t_s^6 + 9*t_s^5)) - (125*t^5*t_s^3*x09)/(2*(500*t_s^6 + 9*t_s^5))];
states = real(states);
end