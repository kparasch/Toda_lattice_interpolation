import numpy as np
import matplotlib.pyplot as plt
import rich.progress

class Particle:
    def __init__(self, q1, p1, q2, p2, q3, p3):
        self.q0 = 0.
        self.q4 = 0.
        self.q1 = q1
        self.p1 = p1
        self.q2 = q2
        self.p2 = p2
        self.q3 = q3
        self.p3 = p3
    
    def x_vec(self):
        return np.array([self.q1, self.p1, self.q2, self.p2, self.q3, self.p3])

    def f_vec(self):
        return np.array([0.5*np.exp(-0.5*(self.q2 - self.q1)), -0.5*self.p1,
                         0.5*np.exp(-0.5*(self.q3 - self.q2)), -0.5*self.p2,
                         0.5*np.exp(-0.5*(self.q4 - self.q3)), -0.5*self.p3])

def Vr(r):
    return np.exp(-r) + r - 1

def Toda_V_3D(q1,q2,q3):
    return Vr(q1) + Vr(q2-q1) + Vr(q3-q2) + V(-q3)

def kinetic_kick(p, dt):
    p.q1 += p.p1 * dt
    p.q2 += p.p2 * dt
    p.q3 += p.p3 * dt
    return

def Toda_V_kick(p, dt):
    p.p1 += dt * (- np.exp( ( p.q1 - p.q2 ) ) )
    #p.p1 += dt * ( np.exp( - ( p.q1 - p.q0 ) ) - np.exp( - ( p.q2 - p.q1 ) ) )
    p.p2 += dt * ( np.exp( - ( p.q2 - p.q1 ) ) - np.exp( - ( p.q3 - p.q2 ) ) )
    p.p3 += dt * ( np.exp( - ( p.q3 - p.q2 ) ) )
    #p.p3 += dt * ( np.exp( - ( p.q3 - p.q2 ) ) - np.exp( - ( p.q4 - p.q3 ) ) )
    return

part = Particle(0, 0., 0., 0., 0., 0.)
turns=100000
dt = 1.e-4
dt2 = dt/2.
TT = np.linspace(0,turns*dt, turns+1)
x_tbt = np.empty([turns+1,6])
f_tbt = np.empty([turns+1,6])

x_tbt[0,:] = part.x_vec()
f_tbt[0,:] = part.f_vec()
for ii in rich.progress.track(range(turns)):
    kinetic_kick(part,dt2)
    Toda_V_kick(part,dt)
    kinetic_kick(part,dt2)
    x_tbt[ii+1,:] = part.x_vec()
    f_tbt[ii+1,:] = part.f_vec()

plt.close('all')
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(321)
ax1.plot(TT, x_tbt[:,0])
ax2 = fig1.add_subplot(322)
ax2.plot(TT, x_tbt[:,1])
ax3 = fig1.add_subplot(323)
ax3.plot(TT, x_tbt[:,2])
ax4 = fig1.add_subplot(324)
ax4.plot(TT, x_tbt[:,3])
ax5 = fig1.add_subplot(325)
ax5.plot(TT, x_tbt[:,4])
ax6 = fig1.add_subplot(326)
ax6.plot(TT, x_tbt[:,5])

fig21 = plt.figure(2)
ax21 = fig21.add_subplot(311)
ax21.plot(x_tbt[:,0], x_tbt[:,1])
ax22 = fig21.add_subplot(312)
ax22.plot(x_tbt[:,2], x_tbt[:,3])
ax23 = fig21.add_subplot(313)
ax23.plot(x_tbt[:,4], x_tbt[:,5])

fig3 = plt.figure(3)
ax31 = fig3.add_subplot(321)
ax31.plot(TT, f_tbt[:,0])
ax32 = fig3.add_subplot(322)
ax32.plot(TT, f_tbt[:,1])
ax33 = fig3.add_subplot(323)
ax33.plot(TT, f_tbt[:,2])
ax34 = fig3.add_subplot(324)
ax34.plot(TT, f_tbt[:,3])
ax35 = fig3.add_subplot(325)
ax35.plot(TT, f_tbt[:,4])
ax36 = fig3.add_subplot(326)
ax36.plot(TT, f_tbt[:,5])

fig41 = plt.figure(4)
ax41 = fig41.add_subplot(311)
ax41.plot(f_tbt[:,0], f_tbt[:,1])
ax42 = fig41.add_subplot(312)
ax42.plot(f_tbt[:,2], f_tbt[:,3])
ax43 = fig41.add_subplot(313)
ax43.plot(f_tbt[:,4], f_tbt[:,5])
plt.show()
