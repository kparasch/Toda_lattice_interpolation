from TricubicInterpolation.cTricubic import Tricubic_Interpolation as TCI
import numpy as np
import matplotlib.pyplot as plt
import sympy

def kick1(q1,q2):
    return (- np.exp( ( q1 - q2 ) ) )


class TricubTodaKick:
    def __init__(self, Nx, Ny, x0=2, y0=2):
        x,y,z = sympy.symbols('x y z')
        f = sympy.exp(x-y)
        dfdx = sympy.diff(f,x)
        dfdy = sympy.diff(f,y)
        dfdz = sympy.diff(f,z)
        dfdxdy = sympy.diff(dfdx,y)
        dfdxdz = sympy.diff(dfdx,z)
        dfdydz = sympy.diff(dfdy,z)
        dfdxdydz = sympy.diff(dfdxdy,z)
        lamf = sympy.lambdify((x,y,z), f, modules='numpy')
        lamdfdx = sympy.lambdify((x,y,z), dfdx, modules='numpy')
        lamdfdy = sympy.lambdify((x,y,z), dfdy, modules='numpy')
        lamdfdz = sympy.lambdify((x,y,z), dfdz, modules='numpy')
        lamdfdxdy = sympy.lambdify((x,y,z), dfdxdy, modules='numpy')
        lamdfdxdz = sympy.lambdify((x,y,z), dfdxdz, modules='numpy')
        lamdfdydz = sympy.lambdify((x,y,z), dfdydz, modules='numpy')
        lamdfdxdydz = sympy.lambdify((x,y,z), dfdxdydz, modules='numpy')

        self.x0 = -abs(x0)
        self.y0 = -abs(x0)
        self.z0 = -50

        self.dx = (2*abs(x0))/(Nx-1)
        self.dy = (2*abs(y0))/(Ny-1)
        self.dz = 10

        Nz = 11

        A = np.empty([Nx,Ny,Nz,8])
        for i in range(Nx):
            xi = self.dx*i + self.x0
            for j in range(Ny):
                yi = self.dy*j + self.y0
                for k in range(Nz):
                    zi = self.dz*k + self.z0
                    A[i,j,k,0] = lamf(xi, yi, zi) 
                    A[i,j,k,1] = lamdfdx(xi, yi, zi)
                    A[i,j,k,2] = lamdfdy(xi, yi, zi)
                    A[i,j,k,3] = lamdfdz(xi, yi, zi)
                    A[i,j,k,4] = lamdfdxdy(xi, yi, zi)
                    A[i,j,k,5] = lamdfdxdz(xi, yi, zi)
                    A[i,j,k,6] = lamdfdydz(xi, yi, zi)
                    A[i,j,k,7] = lamdfdxdydz(xi, yi, zi)


        
        self.TI = TCI(A, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, method='Exact')
        
    def kick_p1(self, q1, q2):
        if q1 > 10 and q2 > 10:
            k = int(q1//10)*10
            return -self.TI.ddx(q1-k, q2-k, 0)
        else:
            return -self.TI.ddx(q1, q2, 0)
    def kick_p2(self, q1, q2):
        if q1 > 10 and q2 > 10:
            k = int(q1//10)*10
            return -self.TI.ddy(q1-k, q2-k, 0)
        else:
            return -self.TI.ddy(q1, q2, 0)
        # return -self.TI.ddy(q1, q2, 0)

# yo=0.
# LTK = TricubTodaKick(11,11)
# xp = np.linspace(-1.5,1.5,1001)
# xpt = np.linspace(-1.5,1.5,1001)
# yp = np.array([LTK.kick_p1(xx, yo) for xx in xp])
# ypt = kick1(xpt, yo)
# 
# plt.plot(xp,yp)
# plt.plot(xpt,ypt)
# plt.show()