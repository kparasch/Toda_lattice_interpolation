from TricubicInterpolation.pyTricubic import Trilinear_Interpolation as TI
import numpy as np
import matplotlib.pyplot as plt

def kick1(q1,q2):
    return (- np.exp( ( q1 - q2 ) ) )


class LinTodaKick:
    def __init__(self, Nx, Ny, x0=2, y0=2):
        self.x0 = -abs(x0)
        self.y0 = -abs(x0)
        self.z0 = -50

        self.dx = (2*abs(x0))/(Nx-1)
        self.dy = (2*abs(y0))/(Ny-1)
        self.dz = 10

        Nz = 11
        A = np.empty([Nx, Ny, Nz])

        for i in range(Nx):
            xi=self.dx*i + x0
            for j in range(Ny):
                yi=self.dy*j + y0
                A[i,j,:] = kick1(xi, yi) 
        
        self.TI = TI(A, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz)
        
    def kick_p1(self, q1, q2):
        if q1 > 10 and q2 > 10:
            k = int(q1//10)*10
            return self.TI.interp(q1-k, q2-k, 0, 0)
        else:
            return self.TI.interp(q1, q2, 0, 0)
    def kick_p2(self, q1, q2):
        if q1 > 10 and q2 > 10:
            k = int(q1//10)*10
            return -self.TI.interp(q1-k, q2-k, 0, 0)
        else:
            return -self.TI.interp(q1, q2, 0, 0)
        # return -self.TI.interp(q1, q2, 0, 0)

# yo=0.
# LTK = LinTodaKick(11,11)
# xp = np.linspace(-1.5,1.5,1001)
# xpt = np.linspace(-1.5,1.5,1001)
# yp = np.array([LTK.TI.interp(xx, yo, 0, 0) for xx in xp])
# ypt = kick1(xpt, yo)
# 
# plt.plot(xp,yp)
# plt.plot(xpt,ypt)
# plt.show()