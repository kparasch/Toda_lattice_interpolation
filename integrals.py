import matplotlib.pyplot as plt
plt.style.use('kostas2')

NN_lin = [11, 31, 51, 101, 201, 501, 1001, 2001, 5001]
NN_cub = [11, 31, 51, 101, 201, 501, 1001, 2001]

H = 1.28e-11
H1 = 1.20e-14
J1 = 1.30e-11
I1 = 1.18e-11

H_lin = [2.49e-2, 2.17e-3, 8.84e-4, 2.20e-4, 5.55e-5, 8.48e-6, 2.15e-6, 5.44e-7, 8.67e-8]
H1_lin = [1.80e-14, 1.40e-14, 1.01e-14, 1.37e-14, 1.78e-14, 8.15e-15, 3.43e-14, 1.17e-14, 3.04e-14]
J1_lin = [9.96e-2, 8.68e-3, 3.53e-3, 8.80e-4, 2.22e-4, 3.39e-5, 8.63e-6,  2.18e-6, 3.47e-7]
I1_lin = [2.51e-2, 2.25e-3, 9.68e-4, 2.46e-4, 6.25e-5, 9.18e-6, 2.37e-6, 6.03e-7, 9.57e-8]

H_cub = [8.67e-5, 1.18e-6, 1.34e-7, 9.76e-9, 5.91e-10, 5.64e-11, 1.39e-11, 1.29e-11]
H1_cub = [8.36e-5, 1.25e-6, 3.85e-7, 4.20e-8, 1.11e-9, 6.64e-11, 1.53e-11, 5.76e-13]
J1_cub = [3.11e-4, 3.40e-6, 1.83e-6, 1.48e-7, 3.79e-9, 3.17e-10, 2.42e-10, 1.29e-11]
I1_cub = [2.76e-4, 4.63e-6, 1.00e-6, 1.17e-7, 3.34e-9, 1.69e-10, 4.09e-11, 1.37e-11]

fig = plt.figure(1, figsize=[16,9])
ax1 = fig.add_subplot(221)
ax1.plot(NN_lin, H_lin, 'ro')
ax1.plot(NN_cub, H_cub, 'bo')
ax1.axhline(H, c='k')
ax1.set_ylabel('std($H$)')

ax2 = fig.add_subplot(222)
ax2.plot(NN_lin, H1_lin, 'ro')
ax2.plot(NN_cub, H1_cub, 'bo')
ax2.axhline(H1, c='k')
ax2.set_ylabel('std($H_1$)')

ax3 = fig.add_subplot(223)
ax3.plot(NN_lin, J1_lin, 'ro')
ax3.plot(NN_cub, J1_cub, 'bo')
ax3.axhline(J1, c='k')
ax3.set_ylabel('std($J_1$)')

ax4 = fig.add_subplot(224)
ax4.plot(NN_lin, I1_lin, 'ro')
ax4.plot(NN_cub, I1_cub, 'bo')
ax4.axhline(I1, c='k')
ax4.set_ylabel('std($I_1$)')


for ax in [ax1, ax2, ax3, ax4]:
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('$N$')
    ax.grid(0)

plt.tight_layout()
plt.show()