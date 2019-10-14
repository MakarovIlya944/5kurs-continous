from Element import Element
from itertools import accumulate, combinations_with_replacement, permutations
from math import floor
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
import operator
from mpl_toolkits.mplot3d import Axes3D
from time import gmtime, strftime
import numpy as np
import logging

with open(f'log-{strftime("%H-%M-%S", gmtime())}.txt','w') as f:
    pass
logging.basicConfig(filename=f'log-{strftime("%H-%M-%S", gmtime())}.txt', level=logging.DEBUG)
logger = logging.getLogger('Main')

numberLocalMatrix = 0

class Spline():

    paint_K = 10

    K = 1
    dim = 1

    w = []

    net = []
    elements = []
    points = []
    h = []
    mn = []
    mx = []
    indexs = []
    n = []

    local_matrix = np.array([ 
        np.array([12, 6, -12, 6]),
        np.array([6, 4, -6, 2]),
        np.array([-12, -6, 12, -6]),
        np.array([6, 2, -6, 4])])

    A = []
    F = []

    __f = []
    __localMatrixes = []

    def __f1(x, l, h):
        t = (x - l) / h
        return 1 - 3*t*t + 2*t*t*t

    def __f2(x, l, h):
        t = (x - l) / h
        return h*(t - 2*t*t + t*t*t)

    def __f3(x, l, h):
        t = (x - l) / h
        return 3*t*t - 2*t*t*t

    def __f4(x, l, h):
        t = (x - l) / h
        return h*(-t*t+t*t*t)

    def __psi(self, el, x, i):
        s = 1
        for k in range(self.dim):
            _i = self.indexs[i][k]
            s *= self.__f[_i](x[k], el.mn[k], self.h[k])
        return s

    def __localMatrix(h):
        local_matrix = Spline.local_matrix * (1/h)

        local_matrix[0][1] /= h
        local_matrix[0][3] /= h
        local_matrix[1][0] /= h
        local_matrix[1][2] /= h
        local_matrix[2][1] /= h
        local_matrix[2][3] /= h
        local_matrix[3][0] /= h
        local_matrix[3][2] /= h

        local_matrix[0][0] /= (h*h)
        local_matrix[0][2] /= (h*h)
        local_matrix[2][0] /= (h*h)
        local_matrix[2][2] /= (h*h)
        return local_matrix

    def __elem(self, ind):
        n = self.n
        p = combinations_with_replacement([0,1],self.dim)
        res = []
        tmp = set()
        for comb in p:
            comb = permutations(comb, self.dim)
            for el in comb:
                tmp.add(el)
        for el in tmp:
            res.append(int(np.sum((np.array(ind) + np.array(el)) * n)))
        res.sort()
        return res

    def __elemInit(self, d, mas):
        if d == 0:
            res = []
            for i in range(self.K):
                mn = list(np.append(mas, i)*self.h)
                mn = np.array(mn)
                res.append(Element(mn + self.mn, self.__elem(np.append(mas, i))))
            return res
        return [self.__elemInit(d-1, np.append(mas, i)) for i in range(self.K)]

    def __elemAdd(self, d, mas):
        if d == 0:
            self.AppendLocalMatrix(mas)
            return
        for i in mas:
            self.__elemAdd(d-1, i)

    def __init__(self, file, _K, _paint_K=10):
        logger.info('Init')

        self.__f = [Spline.__f1, Spline.__f2, Spline.__f3, Spline.__f4]
        self.kElem = _K
        self.paint_K = _paint_K

        points = np.loadtxt(file)
        logger.info(f'Read {points}')

        f = [p[-1] for p in points]
        self.f = f
        a = len(points[0]) - 1
        points = [p[:a] for p in points]
        self.points = points

        mx = []
        dim = len(points[0])
        self.dim = dim
        
        a = [max(points, key=lambda x: x[el])[el] for el in range(dim)]
        mx = np.array(a)
        a = [min(points, key=lambda x: x[el])[el] for el in range(dim)]
        mn = np.array(a)
        self.mn = mn
        self.mx = mx
        K = self.kElem
        self.h = (mx - mn) * (1.0 / K)
        h = self.h
        k_part = 1.01
        mx *= k_part

        self.w = np.ones(len(points))
        
        for _h in self.h:
            self.__localMatrixes.append(Spline.__localMatrix(_h))

        self.kElem = K
        self.kNode = [k+1 for k in K]
        kn = self.kNode
        self.nElem = list(accumulate(K, operator.mul))[-1]
        self.nNodes = list(accumulate([el+1 for el in K], operator.mul))[-1]

        N = list(accumulate(K[:len(K)-1], operator.mul))
        N.insert(0, 1)
        N.reverse()
        self.n = np.array(N)

        if dim == 1:
            self.netX = np.ones(kn[0])*mn + np.array(range(kn[0]))*self.h
            self.net = self.netX
            self.elements = [Element(self.net[i], [i, i+1]) for i in range(self.nElem)]
        elif dim == 2:
            self.netX = np.ones(kn[0])*mn[0] + np.array(range(kn[0]))*self.h[0]
            self.netY = np.ones(kn[1])*mn[1] + np.array(range(kn[1]))*self.h[1]
            mesh = np.meshgrid(self.netX, self.netY)
            f_i = lambda ind: [ind, ind + 1, ind + kn[0], ind + kn[0] + 1]
            for i in range(K[1]):
                for j in range(K[0]):
                    ind = i*kn[0] + j
                    self.net.append(np.array([mesh[0][i][j],mesh[1][i][j]]))
                    self.elements.append(Element(self.net[ind], f_i(ind)))
                self.net.append(np.array([mesh[0][i][K[0]],mesh[1][i][K[0]]]))
            for j in range(len(self.netX)):
                self.net.append(np.array([mesh[0][K[1]][j],mesh[1][K[1]][j]]))

        logger.info(f'{self.nElem} elements created')
        l = list(accumulate([k*2 for k in kn], operator.mul))[-1]
        l = range(l)
        self.A = [np.zeros(len(l)) for i in l] 
        self._A = np.copy(self.A)
        self.F = np.zeros(len(l))

        for I,el in enumerate(points):
            p = [floor((el[i] - mn[i]) / h[i]) for i in range(dim)]
            p = np.array([i if K[v] != i else K[v] - 1 for v,i in enumerate(p)]) * np.array(self.n)
            t = self.elements[np.sum(p)]
            logger.info(f'Point {el} added to element {t.i}')
            t.addP(el)
            t.addF(f[I])
            t.addW(self.w[I])

        with open(f'{dim}d.txt','r') as f:
            lines = f.readlines()
            # Not implemented dim>10
            self.indexs = [[int(c)-1 for c in str(int(l))] for l in lines]

    def Calculate(self):
        logger.info('Calculate')
        self.MakeMatrix()
        self.Solve()

    def MakeMatrix(self):
        logger.info('MakeMatrix')
        for el in self.elements:
            self.AppendLocalMatrix(el)

    def AppendLocalMatrix(self, el):
        global numberLocalMatrix
        logger.info('MakeLocalMatrix')
        dim = self.dim
        nums = range(pow(4, dim))
        psi = self.__psi
        local_f_number = 2**dim

        L = self.__localMatrixes

        inds = self.indexs
        for I in nums:
            for J in nums:
                value = el.b
                for i in range(self.dim):
                    logger.debug(f'L[{i}][{inds[I][i]}][{inds[J][i]}] = {L[i][inds[I][i]][inds[J][i]]}')
                    value *= L[i][inds[I][i]][inds[J][i]]

                for i, p in enumerate(el.p):
                    value += el.w[i] * psi(el, p, I) * psi(el, p, J)
                    # print(f'el:{numberLocalMatrix}n:{i}value:{value}')
                
                i = el.nodes[I // local_f_number] * local_f_number + I % local_f_number
                j = el.nodes[J // local_f_number] * local_f_number + J % local_f_number
                logger.debug(f'i={i}\tj={j}')
                
                self._A[i][j] = numberLocalMatrix + 1
                # self._A[i][j] = value
                self.A[i][j] += value

            for _i, p in enumerate(el.p):
                self.F[i] += el.w[_i] * psi(el, p, I) * el.f[_i]
        np.savetxt(f'step_{numberLocalMatrix}.txt',self._A,fmt='%.0f')
        numberLocalMatrix += 1

    def Solve(self):
        logger.info('Solve')
        self.answer = np.linalg.solve(self.A, self.F)
        return self.answer

    def Paint(self):
        logger.info('Paint')
        x = []
        y = []
        z = []
        K = self.paint_K
        psi = self.__psi
        borders = []

        # local functions per node count
        lfnn = 2**self.dim

        # local functions per element count 
        lfne = 4**self.dim

        # range local functions per element count
        rle = range(lfne)

        if self.dim == 1:
            elem_steps = [self.h[0] * (el/K) for el in range(K)]
            for i,el in enumerate(self.elements):
                _x = [_el + el.mn for _el in elem_steps]
                x.extend(_x)
                y.extend([list(accumulate([self.answer[v+i*lfnn] * psi(el, y, v) for v in rle]))[-1] for y in _x])
                borders.append(el.mn)
            _x = el.mn + self.h[0]
            x.append(_x)
            y.append(list(accumulate([self.answer[v+i*lfnn] * psi(el, _x, v) for v in rle]))[-1])
            plt.plot(x, y, '-', self.points, self.f, 'o', borders, np.ones(self.K) * 5, '+')
            plt.show()
        elif self.dim == 2:
            fig = plt.figure()
            ax = fig.gca(projection='3d')   
            elem_steps = []
            z = [np.zeros(len(self.elements) * K) for el in range(len(self.elements) * K)]
            for i in range(self.dim):
                elem_steps.append(list(accumulate(np.ones(K) * (self.h[i] / K))))
            for i,el in enumerate(self.elements):
                _x = [_el + el_x[0].mn[0] for _el in elem_steps[0]]
                x.extend(_x)
                for j,el_y in enumerate(el_x):
                    _y = [_el + el_y.mn[1] for _el in elem_steps[1]]
                    if i == 0:
                        y.extend(_y)
                    for cur_x in range(K):
                        for cur_y in range(K):
                            _I = i*K + cur_x
                            _J = j*K + cur_y
                            z[_I][_J] = np.sum([self.answer[el_y.nodes[v//lfnn]*lfnn+v%lfnn] * psi(el_y, [x[_I],y[_J]], v) for v in rle])
            x, y = np.meshgrid(x, y)
            z = np.array(z)
            surf = ax.plot_surface(x, y, z)
            # ax.hold(True)
            # ax = fig.add_subplot(111, projection='3d')
            # for i,p in enumerate(self.points):
            #     # Make data
            #     u = np.linspace(0, 2 * np.pi, 100)
            #     v = np.linspace(0, np.pi, 100)
            #     x = 10 * np.outer(np.cos(u), np.sin(v))
            #     y = 10 * np.outer(np.sin(u), np.sin(v))
            #     z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))
            #     ax.plot_surface(x, y, z, color='b')
            # a = [np.append(p, self.f[i]) for i,p in enumerate(self.points)]
            # ax.scatter(a, color='green')
            plt.show()
