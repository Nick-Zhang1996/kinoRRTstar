# DI3d_* function in matlab
# generated to solve optimal control problems
import numpy as np

class DoubleIntegrator:
    # the naming of parameters follow the original implementation
    # x0? dentoes the initial state
    # x1? denotes the final state
    # if velocity at final state is free, there are only 3 states
    # otherwise x,y,z,vx,vy,vz are specified for each state, in that order
    def __init__(self):
        return
    
    # following functions return a single scalar cost
    # cost from state a to b
    def DI3d_costFreeVel(self,t_s,x01,x02,x03,x04,x05,x06,x11,x12,x13):
        return t_s + (3*x04**2 + 3*x05**2 + 3*x06**2)/t_s + (6*x01*x04 + 6*x02*x05 + 6*x03*x06 - 6*x04*x11 \
        - 6*x05*x12 - 6*x06*x13)/t_s**2 + (3*x01**2 - 6*x01*x11 + 3*x02**2 - 6*x02*x12 + 3*x03**2 - 6*x03*x13 + 3*x11**2 + 3*x12**2 + 3*x13**2)/t_s**3

    def DI3d_cost(self,t_s,x01,x02,x03,x04,x05,x06,x11,x12,x13,x14,x15,x16):
        return t_s + (12*x01**2)/t_s**3 + (4*x04**2)/t_s + (12*x02**2)/t_s**3 + (4*x05**2)/t_s + (12*x03**2)/t_s**3 + (4*x06**2)/t_s \
               + (12*x11**2)/t_s**3 + (4*x14**2)/t_s + (12*x12**2)/t_s**3 + (4*x15**2)/t_s + (12*x13**2)/t_s**3 + (4*x16**2)/t_s  \
               + (12*x01*x04)/t_s**2 + (12*x02*x05)/t_s**2 + (12*x03*x06)/t_s**2 - (24*x01*x11)/t_s**3 + (12*x01*x14)/t_s**2 \
               - (12*x04*x11)/t_s**2 - (24*x02*x12)/t_s**3 + (4*x04*x14)/t_s + (12*x02*x15)/t_s**2 - (12*x05*x12)/t_s**2 \
               - (24*x03*x13)/t_s**3 + (4*x05*x15)/t_s + (12*x03*x16)/t_s**2 - (12*x06*x13)/t_s**2 + (4*x06*x16)/t_s - (12*x11*x14)/t_s**2 - (12*x12*x15)/t_s**2 - (12*x13*x16)/t_s**2

    # following functions calculate state vector from initial to final state
    # n_interior : interior points, return vector has interior+2 entries
    # return vector dim: (interior+2, 6)
    # NOTE this is a little different from original matlab implementation in that
    # it returns the initial and final state as well
    # it is not necessary to check boundary states for collision (should already be checked previously)
    # however, this allows user to get velocity at final state
    def DI3d_stateFreeVel(self,t_s,x01,x02,x03,x04,x05,x06,x11,x12,x13,n_interior=10):
        t = np.linspace(0,t_s,n_interior+2)
        p = np.vstack([  x01 + t*x04 + (t**2*(t - 3*t_s)*(x01 - x11 + t_s*x04))/(2*t_s**3),
                    x02 + t*x05 + (t**2*(t - 3*t_s)*(x02 - x12 + t_s*x05))/(2*t_s**3),
                    x03 + t*x06 + (t**2*(t - 3*t_s)*(x03 - x13 + t_s*x06))/(2*t_s**3),
                    x04 + (3*t*(t - 2*t_s)*(x01 - x11 + t_s*x04))/(2*t_s**3),
                    x05 + (3*t*(t - 2*t_s)*(x02 - x12 + t_s*x05))/(2*t_s**3),
                    x06 + (3*t*(t - 2*t_s)*(x03 - x13 + t_s*x06))/(2*t_s**3)])

        return p.T

    def DI3d_state(self,t_s,x01,x02,x03,x04,x05,x06,x11,x12,x13,x14,x15,x16,n_interior=10):
        t = np.linspace(0,t_s,n_interior+2)
        p = np.vstack([x11 + x14*(t - t_s) + ((t - t_s)**3*(2*x01 - 2*x11 + t_s*x04 + t_s*x14))/t_s**3 + ((t - t_s)**2*(3*x01 - 3*x11 + t_s*x04 + 2*t_s*x14))/t_s**2,
              x12 + x15*(t - t_s) + ((t - t_s)**3*(2*x02 - 2*x12 + t_s*x05 + t_s*x15))/t_s**3 + ((t - t_s)**2*(3*x02 - 3*x12 + t_s*x05 + 2*t_s*x15))/t_s**2,
              x13 + x16*(t - t_s) + ((t - t_s)**3*(2*x03 - 2*x13 + t_s*x06 + t_s*x16))/t_s**3 + ((t - t_s)**2*(3*x03 - 3*x13 + t_s*x06 + 2*t_s*x16))/t_s**2,
              x14 + (2*(t - t_s)*(3*x01 - 3*x11 + t_s*x04 + 2*t_s*x14))/t_s**2 + (3*(t - t_s)**2*(2*x01 - 2*x11 + t_s*x04 + t_s*x14))/t_s**3,
              x15 + (2*(t - t_s)*(3*x02 - 3*x12 + t_s*x05 + 2*t_s*x15))/t_s**2 + (3*(t - t_s)**2*(2*x02 - 2*x12 + t_s*x05 + t_s*x15))/t_s**3,
              x16 + (2*(t - t_s)*(3*x03 - 3*x13 + t_s*x06 + 2*t_s*x16))/t_s**2 + (3*(t - t_s)**2*(2*x03 - 2*x13 + t_s*x06 + t_s*x16))/t_s**3])

        return p.T

    # following functions solve the particular polynomial equation and
    # return the smallest positive real root
    def DI3d_timeFreeVel(self,x01, x02, x03, x04, x05, x06, x11, x12, x13):
        p = [ 1, 0, 
            - 3*x04**2 - 3*x05**2 - 3*x06**2, 12*x04*x11 - 12*x02*x05 - 12*x03*x06 - 12*x01*x04 + 12*x05*x12 + 12*x06*x13, 
            - 9*x01**2 + 18*x01*x11 - 9*x02**2 + 18*x02*x12 - 9*x03**2 + 18*x03*x13 - 9*x11**2 - 9*x12**2 - 9*x13**2]
        x = np.roots(p)
        mask = np.logical_and(x.imag == 0, x.real > 0)
        val = np.min(x[mask].real)
        #if (np.abs(self.quartic_root(p) - val) > 0.01):
        #    breakpoint()
        return val

    def DI3d_time(self,x01, x02, x03, x04, x05, x06, x11, x12, x13, x14, x15, x16):
        p = [ 1, 0, 
        - 4*x04**2 - 4*x04*x14 - 4*x05**2 - 4*x05*x15 - 4*x06**2 - 4*x06*x16 - 4*x14**2 - 4*x15**2 - 4*x16**2, 
        24*x04*x11 - 24*x02*x05 - 24*x03*x06 - 24*x01*x14 - 24*x01*x04 - 24*x02*x15 + 24*x05*x12 - 24*x03*x16 + 24*x06*x13 + 24*x11*x14 + 24*x12*x15 + 24*x13*x16, 
        - 36*x01**2 + 72*x01*x11 - 36*x02**2 + 72*x02*x12 - 36*x03**2 + 72*x03*x13 - 36*x11**2 - 36*x12**2 - 36*x13**2]

        x = np.roots(p)
        mask = np.logical_and(x.imag == 0, x.real > 0)
        val = np.min(x[mask].real)
        #if (np.abs(self.quartic_root(p) - val) > 0.01):
        #    breakpoint()
        return val

    # return the smallest real root of quartic polynomial
    # according to wikipedia
    def quartic_root(self,coeffs):
        a,b,c,d,e = coeffs
        D0 = c*c - 3*b*d + 12*a*e
        D1 = 2*c*c*c - 9*b*c*d + 27*b*b*e + 27*a*d*d - 72*a*c*e
        Q = ((D1 + (D1*D1 - 4*D0*D0*D0)**0.5)/2)**(1.0/3.0)
        p = (8*a*c - 3*b*b)/(8*a*a)
        q = (b*b*b - 4*a*b*c + 8*a*a*d)/(8*a*a*a)
        S = 0.5 * (-2.0/3.0 * p + 1.0/(3*a)*(Q + D0/Q))**0.5

        val = None

        temp1 = -4*S*S - 2*p + q/S
        if (temp1 > 0):
            candidate = -b/(4*a) - S + 0.5*temp1**0.5
            if (candidate > 0):
                if (val is None):
                    val = candidate
                elif (candidate < val):
                    val = candidate

            candidate = -b/(4*a) - S - 0.5*temp1**0.5
            if (candidate > 0):
                if (val is None):
                    val = candidate
                elif (candidate < val):
                    val = candidate

        temp2 = -4*S*S - 2*p - q/S
        if (temp2 > 0):
            candidate = -b/(4*a) + S + 0.5*temp2**0.5
            if (candidate > 0):
                if (val is None):
                    val = candidate
                elif (candidate < val):
                    val = candidate

            candidate = -b/(4*a) + S - 0.5*temp2**0.5
            if (candidate > 0):
                if (val is None):
                    val = candidate
                elif (candidate < val):
                    val = candidate
        if (val is None):
            breakpoint()
        return val


if __name__=="__main__":
    DI = DoubleIntegrator()

    # should be 5.26
    ans = DI.DI3d_time(17.92223,5.97527,8.84017,-0.60555,-0.61722,1.36524,18.00000,8.00000,8.00000,0.00000,0.00000,0.00000)
    print(ans-5.26)

    # should be 4.6188
    cost_free = DI.DI3d_costFreeVel(3.46410,2.00000,2.00000,2.00000,0.00000,0.00000,0.00000,5.24243,3.90719,3.35989)
    print(cost_free-4.6188)
    tf = DI.DI3d_timeFreeVel(2.00000,2.00000,2.00000,0.00000,0.00000,0.00000,5.24243,3.90719,3.35989)
    print(tf-3.46410)
    # now get cost
    state = DI.DI3d_stateFreeVel(tf,2.00000,2.00000,2.00000,0.00000,0.00000,0.00000,5.24243,3.90719,3.35989)
    final_state = state[-1,:]
    cost_fixed = DI.DI3d_cost(tf,2.00000,2.00000,2.00000,0.00000,0.00000,0.0,*final_state)
    print(cost_fixed-cost_free)


    # should be
    expected = \
    [2.0338,    2.1312,    2.2856,    2.4908,    2.7405,    3.0282,    3.3476,    3.6925,    4.0564,    4.4329,
    2.0231,    2.0894,    2.1947,    2.3347,    2.5049,    2.7011,    2.9189,    3.1540,    3.4021,    3.6589,
    2.0252,    2.0975,    2.2124,    2.3650,    2.5506,    2.7646,    3.0021,    3.2586,    3.5292,    3.8092,
    0.2116,    0.4031,    0.5744,    0.7255,    0.8565,    0.9674,    1.0581,    1.1286,    1.1790,    1.2092,
    0.1443,    0.2748,    0.3916,    0.4947,    0.5840,    0.6596,    0.7214,    0.7695,    0.8039,    0.8245,
    0.1574,    0.2997,    0.4271,    0.5395,    0.6369,    0.7194,    0.7868,    0.8393,    0.8767,    0.8992]
    expected = np.array(expected).reshape(6,10)
    

    ans = DI.DI3d_stateFreeVel(3.46410,2.00000,2.00000,2.00000,0.00000,0.00000,0.00000,4.81585,3.91998,4.09397,)
    error = (ans[1:-1]-expected.T)
    print(np.linalg.norm(error))
