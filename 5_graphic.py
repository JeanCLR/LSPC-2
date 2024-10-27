#!/usr/bin/python3
##################################################################################
import numpy as np
import matplotlib.pyplot as plt
import math
##################################################################################
import warnings
warnings.filterwarnings("ignore")
##################################################################################
##################################################################################

radio = np.loadtxt('file1.dat',float)
theta = np.loadtxt('file2.dat',float)
temp = np.loadtxt('file4.dat',float)

fig = plt.figure()
#fig.patch.set_facecolor('red')
ax = plt.subplot(111,projection='polar')

th, r = np.meshgrid(theta,radio)

im = ax.pcolormesh(th,r,temp,cmap='jet',shading='nearest')
cbar = plt.colorbar(im, aspect=15, pad=0.18)#pad=0.17, fraction=0.07, aspect=8)
cbar.ax.set_title(r'$T(r,\theta)\ (^\circ C)$')
cbar.ax.get_yaxis().labelpad = 18
cbar.ax.yaxis.tick_left()
ax.axis('off')
#ax.set_thetagrids([])
#ax.set_rgrids(radio)
#ax.set_rlabel_position(55)
ax.tick_params('both', labelsize=12)
ax.grid(color='black', alpha=0.7, linewidth=0.15)

#plt.savefig('sol_numerica.eps', bbox_inches='tight')
#plt.savefig('sol_numerica.pdf', bbox_inches='tight')
plt.show()