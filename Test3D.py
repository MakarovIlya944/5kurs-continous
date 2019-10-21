import numpy as np
from matplotlib import pyplot as plt
from celluloid import Camera
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')
camera = Camera(fig)

x = np.arange(0,10,0.1)
y = np.arange(0,10,0.1)
t = np.arange(0,10,1)

xv, yv = np.meshgrid(x, y)
for i in t:
    T = np.array([np.ones(len(x)) for i in  range(len(y))]) * i
    z = np.sin(xv + yv + T)
    ax.plot_surface(xv, yv, z, color='blue')
    camera.snap()
    
animation = camera.animate()
animation.save('celluloid_subplots.gif', writer = 'imagemagick')