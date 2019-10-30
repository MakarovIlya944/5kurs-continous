import numpy as np
from matplotlib import pyplot as plt
from celluloid import Camera
from matplotlib import cm
from itertools import accumulate, combinations_with_replacement, permutations
from mpl_toolkits.mplot3d import Axes3D
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Painter')

class Painter():

    paint_K = 0

    elements = []
    clearPoints = []
    noisePoints = []
    h = []
    
    mx = 0

    dim = 0
    psi = 0
    answer = []

    isPipe = False

    def __init__(self, file, dim, pipeline=False, psi=None, k=None, elements=None, mx=None, kElem=None, h=None, noisePoints=None, clearPoints=None):
        self.answer = np.loadtxt(file)
        self.psi = psi
        self.paint_K = k
        self.dim = dim
        self.noisePoints = noisePoints
        self.clearPoints = clearPoints
        self.elements = elements
        self.mx = mx
        self.kElem = kElem
        self.h = h
        self.isPipe = pipeline

        # local functions per node count
        self.lfnn = 2**self.dim

        # local functions per element count 
        self.lfne = 4**self.dim

        # range local functions per element count
        self.rle = range(self.lfne)

    def __paint1D(self, points=False):
        K = self.paint_K
        psi = self.psi
        
        x = []
        f = []
        lfnn = self.lfnn
        rle = self.rle

        elem_steps = [self.h[0] * (el/K) for el in range(K)]
        for i,el in enumerate(self.elements):
            _x = [_el + el.mn for _el in elem_steps]
            x.extend(_x)
            f.extend([list(accumulate([self.answer[v+i*lfnn] * psi(el, y, v) for v in rle]))[-1] for y in _x])
        _x = el.mn + self.h[0]
        x.append(_x)
        f.append(list(accumulate([self.answer[v+i*lfnn] * psi(el, _x, v) for v in rle]))[-1])
        plt.plot(x, f, '-')
        if points:
            xs = []
            ys = []
            for p in self.clearPoints:
                xs.append(p[0])
                ys.append(p[1])
            plt.plot(xs, ys, 'o', color=(1,0,0))
            xs = []
            ys = []
            for p in self.noisePoints:
                xs.append(p[0])
                ys.append(p[1])
            plt.plot(xs, ys, '*', color=(0,0,1))

    def __paint2D(self, points=False):
        K = self.paint_K
        psi = self.psi

        fig = plt.figure()
        ax = fig.gca(projection='3d')   

        elem_steps = []
        x = []
        y = []
        lfnn = self.lfnn
        rle = self.rle
        rng = range(K+1)

        z = [np.zeros(self.kElem[1] * K + 1) for el in range(self.kElem[0] * K + 1)]
        for i in range(self.dim):
            tmp = np.ones(K) * (self.h[i] / K)
            tmp = np.insert(tmp, 0, 0)[:-1]
            tmp = list(accumulate(tmp))
            elem_steps.append(tmp)

        for el_x in self.elements:
            _x = [_el + el_x[0].mn[0] for _el in elem_steps[0]]
            x.extend(_x)
        x.append(self.mx[0])
        for el_y in self.elements[0]:
            _y = [_el + el_y.mn[1] for _el in elem_steps[1]]
            y.extend(_y)
        y.append(self.mx[1])

        for i,el_x in enumerate(self.elements):
            for j,el_y in enumerate(el_x):
                for cur_x in rng:
                    for cur_y in rng:
                        _I = i*K + cur_x
                        _J = j*K + cur_y
                        z[_I][_J] = np.sum([self.answer[el_y.nodes[v//lfnn]*lfnn+v%lfnn] * psi(el_y, [x[_I],y[_J]], v) for v in rle])
        y, x = np.meshgrid(y, x)
        z = np.array(z)
        ax.plot_surface(x, y, z)
        if points:
            xs = []
            ys = []
            zs = []
            for p in self.clearPoints:
                xs.append(p[0])
                ys.append(p[1])
                zs.append(p[2])
            ax.scatter(xs, ys, zs, marker='o', color=(1,0,0)) 
            xs = []
            ys = []
            zs = []
            for p in self.noisePoints:
                xs.append(p[0])
                ys.append(p[1])
                zs.append(p[2])
            ax.scatter(xs, ys, zs, marker='*', color=(0,0,1)) 
        
    def __paint3D(self, points=False, fromFile=False, testFunc=None, gif=True):

        fig = plt.figure()
        ax = fig.gca(projection='3d')   
        camera = Camera(fig)

        if not fromFile:
            K = self.paint_K
            psi = self.psi
            elem_steps = []
            x = []
            y = []
            z = []
            f = []
            lfnn = self.lfnn
            rle = self.rle
            rng = range(K+1)

            f = [[np.zeros(self.kElem[1] * K + 1) for el in range(self.kElem[0] * K + 1)] for el in range(self.kElem[2] * K + 1)]
            for i in range(self.dim):
                tmp = np.ones(K) * (self.h[i] / K)
                tmp = np.insert(tmp, 0, 0)[:-1]
                tmp = list(accumulate(tmp))
                elem_steps.append(tmp)

            for el_x in self.elements:
                _x = [_el + el_x[0][0].mn[0] for _el in elem_steps[0]]
                x.extend(_x)
            x.append(self.mx[0])
            for el_y in self.elements[0]:
                _y = [_el + el_y[0].mn[1] for _el in elem_steps[1]]
                y.extend(_y)
            y.append(self.mx[1])
            for el_z in self.elements[0][0]:
                _z = [_el + el_z.mn[2] for _el in elem_steps[1]]
                z.extend(_z)
            z.append(self.mx[2])

            for i,el_x in enumerate(self.elements):
                for j,el_y in enumerate(el_x):
                    for k,el_z in enumerate(el_y):
                        for cur_x in rng:
                            _I = i*K + cur_x
                            for cur_y in rng:
                                _J = j*K + cur_y
                                for cur_z in rng:
                                    _K = k*K + cur_z
                                    if testFunc:
                                        f[_K][_I][_J] = testFunc(x[_I],y[_J],z[_K])
                                    else:   
                                        f[_K][_I][_J] = np.sum([self.answer[el_z.nodes[v//lfnn]*lfnn+v%lfnn] * psi(el_z, [x[_I],y[_J],z[_K]], v) for v in rle])
            y, x = np.meshgrid(y, x)
            Z = z
            i = 0
            _____i = 0
            for _f in f:
                z = np.array(_f)
                ax.plot_surface(x, y, z, color='grey')
                if points:
                    xs = []
                    ys = []
                    zs = []
                    for p in self.clearPoints:
                        if abs(p[2] - Z[i]) < 1E-1:
                            xs.append(p[0])
                            ys.append(p[1])
                            zs.append(p[3])
                    ax.scatter(xs, ys, zs, marker='o', color=(1,0,0)) 
                    xs = []
                    ys = []
                    zs = []
                    for p in self.noisePoints:
                        if abs(p[2] - Z[i]) < 1E-1:
                            xs.append(p[0])
                            ys.append(p[1])
                            zs.append(p[3])
                    ax.scatter(xs, ys, zs, marker='*', color=(0,0,1)) 
                if logger.level == logging.DEBUG or self.isPipe:
                    np.savetxt(f'data/z{i}.txt', z, fmt='%1.2f')
                    i += 1
                camera.snap()
                fig.savefig(f'f{i}')
            animation = camera.animate()
            animation.save('data/3d.gif', writer = 'imagemagick')
        else:
            f = []
            for i in range(11):
                f.append(np.loadtxt(f'f{i}.txt'))
            x = np.loadtxt('x.txt')
            y = np.loadtxt('y.txt')
            for i in range(11):
                z = np.array(f[i])
                ax.plot_surface(x, y, z)
                fig.savefig(f'fff{i}')

    def Paint(self, points=False, testFunc=None):
        if self.dim == 1:
            self.__paint1D(points)
            plt.show()
        elif self.dim == 2:
            self.__paint2D(points)
            plt.show()
        elif self.dim == 3:
            self.__paint3D(points,testFunc=testFunc)
        
if __name__ == "__main__":
    p = Painter('answer.txt', 3)
    p.Paint()