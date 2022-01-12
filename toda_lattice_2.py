import numpy as np
import matplotlib.pyplot as plt
import rich.progress
from linear_kick import LinTodaKick
from tricub_kick import TricubTodaKick

class Particle:
    def __init__(self, q1, p1, q2, p2):
        self.q0 = 0.
        self.q3 = 0.
        self.q1 = q1
        self.p1 = p1
        self.q2 = q2
        self.p2 = p2
    
    def x_vec(self):
        return np.array([self.q1, self.p1, self.q2, self.p2])

    def f_vec(self):
        return np.array([0.5*np.exp(-0.5*(self.q2 - self.q1)), -0.5*self.p1,
                         0.5*np.exp(-0.5*(self.q3 - self.q2)), -0.5*self.p2
                         ])

    def I_vec(self):
        H = 0.5*self.p1**2 + 0.5*self.p2**2 + np.exp(self.q1-self.q2) 
        H1 = -0.5*(self.p1 + self.p2)
        J1 = (self.p1 - self.p2)**2 + 4*np.exp(self.q1-self.q2)
        I1 = (self.p1 - self.p2 + np.sqrt(J1))/(self.p1 - self.p2 -np.sqrt(J1)) * np.exp( np.sqrt(J1)*(self.q1 + self.q2)/(self.p1+self.p2) )
        return np.array([H, H1, J1, I1])

def Vr(r):
    return np.exp(-r) + r - 1

def Toda_V_3D(q1,q2,q3):
    return Vr(q1) + Vr(q2-q1) + Vr(q3-q2) + V(-q3)

def kinetic_kick(p, dt):
    p.q1 += p.p1 * dt
    p.q2 += p.p2 * dt
    return

def Toda_V_kick(p, dt):
    p.p1 += dt * (- np.exp( ( p.q1 - p.q2 ) ) )
    p.p2 += dt * (+ np.exp( ( p.q1 - p.q2 ) ) )
    #p.p1 += dt * ( np.exp( - ( p.q1 - p.q0 ) ) - np.exp( - ( p.q2 - p.q1 ) ) )
    #p.p2 += dt * ( np.exp( - ( p.q2 - p.q1 ) ) - np.exp( - ( p.q3 - p.q2 ) ) )
    return

def Toda_Lin_kick(p, dt, LTK):
    p.p1 += dt * LTK.kick_p1(p.q1, p.q2)
    p.p2 += dt * LTK.kick_p2(p.q1, p.q2)
    return

def Toda_Cub_kick(p, dt, CTK):
    p.p1 += dt * CTK.kick_p1(p.q1, p.q2)
    p.p2 += dt * CTK.kick_p2(p.q1, p.q2)
    return


NN = 201
LTK  = LinTodaKick(NN,NN)
#CTK  = LTK
CTK  = TricubTodaKick(NN,NN)

part = Particle(0.5, 1., -0.5, 0.)
part_lin = Particle(0.5, 1., -0.5, 0.)
part_cub = Particle(0.5, 1., -0.5, 0.)

turns=1000000
dt = 1.e-6
dt2 = dt/2.
TT = np.linspace(0,turns*dt, turns+1)
x_tbt = np.empty([turns+1,4])
x_tbt_lin = np.empty([turns+1,4])
x_tbt_cub = np.empty([turns+1,4])
I_tbt = np.empty([turns+1,4])
I_tbt_lin = np.empty([turns+1,4])
I_tbt_cub = np.empty([turns+1,4])

x_tbt[0,:] = part.x_vec()
x_tbt_lin[0,:] = part_lin.x_vec()
x_tbt_cub[0,:] = part_cub.x_vec()
I_tbt[0,:] = part.I_vec()
I_tbt_lin[0,:] = part_lin.I_vec()
I_tbt_cub[0,:] = part_cub.I_vec()
for ii in rich.progress.track(range(turns)):
    kinetic_kick(part,dt2)
    Toda_V_kick(part,dt)
    kinetic_kick(part,dt2)
    x_tbt[ii+1,:] = part.x_vec()
    I_tbt[ii+1,:] = part.I_vec()

    kinetic_kick(part_lin,dt2)
    Toda_Lin_kick(part_lin,dt, LTK)
    kinetic_kick(part_lin,dt2)
    x_tbt_lin[ii+1,:] = part_lin.x_vec()
    I_tbt_lin[ii+1,:] = part_lin.I_vec()

    kinetic_kick(part_cub,dt2)
    Toda_Cub_kick(part_cub,dt, CTK)
    kinetic_kick(part_cub,dt2)
    x_tbt_cub[ii+1,:] = part_cub.x_vec()
    I_tbt_cub[ii+1,:] = part_cub.I_vec()


plt.close('all')
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(221)
ax1.plot(TT, x_tbt[:,0])
ax1.plot(TT, x_tbt_lin[:,0])
ax1.plot(TT, x_tbt_cub[:,0])
ax2 = fig1.add_subplot(222)
ax2.plot(TT, x_tbt[:,1])
ax2.plot(TT, x_tbt_lin[:,1])
ax2.plot(TT, x_tbt_cub[:,1])
ax3 = fig1.add_subplot(223)
ax3.plot(TT, x_tbt[:,2])
ax3.plot(TT, x_tbt_lin[:,2])
ax3.plot(TT, x_tbt_cub[:,2])
ax4 = fig1.add_subplot(224)
ax4.plot(TT, x_tbt[:,3])
ax4.plot(TT, x_tbt_lin[:,3])
ax4.plot(TT, x_tbt_cub[:,3])

fig21 = plt.figure(2)
ax21 = fig21.add_subplot(211)
ax21.plot(x_tbt[:,0], x_tbt[:,1])
ax21.plot(x_tbt_lin[:,0], x_tbt_lin[:,1])
ax21.plot(x_tbt_cub[:,0], x_tbt_cub[:,1])
ax22 = fig21.add_subplot(212)
ax22.plot(x_tbt[:,2], x_tbt[:,3])
ax22.plot(x_tbt_lin[:,2], x_tbt_lin[:,3])
ax22.plot(x_tbt_cub[:,2], x_tbt_cub[:,3])

fig3 = plt.figure(3)
ax31 = fig3.add_subplot(411)
ax31.plot(TT, I_tbt[:,0])
ax31.plot(TT, I_tbt_lin[:,0])
ax31.plot(TT, I_tbt_cub[:,0])
ax32 = fig3.add_subplot(412)
ax32.plot(TT, I_tbt[:,1])
ax32.plot(TT, I_tbt_lin[:,1])
ax32.plot(TT, I_tbt_cub[:,1])
ax33 = fig3.add_subplot(413)
ax33.plot(TT, I_tbt[:,1])
ax33.plot(TT, I_tbt_lin[:,2])
ax33.plot(TT, I_tbt_cub[:,2])
ax34 = fig3.add_subplot(414)
ax34.plot(TT, I_tbt[:,1])
ax34.plot(TT, I_tbt_lin[:,3])
ax34.plot(TT, I_tbt_cub[:,3])

print(f"Exact rms(IoM): {np.std(I_tbt, axis=0)}")
print(f"Linear rms(IoM): {np.std(I_tbt_lin, axis=0)}")
print(f"Cubic rms(IoM): {np.std(I_tbt_cub, axis=0)}")
plt.show()
