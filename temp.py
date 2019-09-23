import numpy as np
import numpy.random
import matplotlib.pyplot as plt


x_max = 10
y_max = 8

dx = 0.1
dy = 0.1

x_edge = np.arange(dx, x_max, dx)
y_edge = np.arange(0, y_max, dy)

# Fr = z / (R - r + z)
def d(r, z):
    if r > z:
        if (z / y_max) > (1 - r / x_max):
            return z / (r * (x_max - r + z)) + r*z/(x_max - r + z)
        else:
            return 1
    return 0

def main():
    print('')

    extent = [0, x_max, 0, y_max]

    heatmap = []
    for x in x_edge:
        tmp = np.array([d(x,y) for y in y_edge])
        heatmap.append(tmp)

    heatmap = np.array(heatmap)

    # Plot heatmap
    plt.clf()
    plt.ylabel('y')
    plt.xlabel('x')
    plt.imshow(heatmap, extent=extent)
    plt.show()


if __name__ == "__main__":
    main()