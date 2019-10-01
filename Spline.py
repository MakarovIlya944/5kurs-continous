from Element import Element
from itertools import accumulate, permutations
from math import floor
import matplotlib.pyplot as plt
import operator
from time import gmtime, strftime
import numpy as np
import logging

with open(f'log-{strftime("%H-%M-%S", gmtime())}.txt','w') as f:
    pass
logging.basicConfig(filename=f'log-{strftime("%H-%M-%S", gmtime())}.txt', level=logging.DEBUG)
logger = logging.getLogger('Main')

class Spline():

    K = 1
    dim = 1

    w = []

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
        p = permutations([0,1],self.dim)
        res = []
        for el in p:
            res.append(np.sum((np.array(ind) + np.array(el)) * n))
        return res

    def __elemInit(self, d, mas):
        if d == 0:
            return [Element(np.append(mas, i)*self.h + self.mn, self.__elem(np.append(mas, i))) for i in range(self.K)]
        return [self.__elemInit(d-1, np.append(mas, i)) for i in range(self.K)]

    def __elemAdd(self, d, mas):
        if d == 0:
            self.AppendLocalMatrix(mas)
            return
        for i in mas:
            self.__elemAdd(d-1, i)

    def __init__(self, file):
        logger.info('Init')

        self.__f = [Spline.__f1, Spline.__f2, Spline.__f3, Spline.__f4]

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
        K = self.K
        self.h = (mx - mn) * (1.0 / K)
        h = self.h
        k_part = 1.01
        mx *= k_part

        self.w = np.ones(len(points))
        
        for _h in self.h:
            self.__localMatrixes.append(Spline.__localMatrix(_h))

        N = list(accumulate([K for e in range(dim-1)], operator.mul))
        N.insert(0, 1)
        self.n = np.array(N)

        self.elements = self.__elemInit(dim-1, [])
        logger.info(f'{K}^{dim} elements created')

        l = range((2*(K+1))**dim)
        self.A = [np.zeros(len(l)) for i in l] 
        self.F = np.zeros(len(l))

        logger.info('-' * 15)
        for I,el in enumerate(points):
            p = [floor((el[i] - mn[i]) / h[i]) for i in range(dim)]
            p = list([i if K != i else K - 1 for i in p])
            t = self.elements
            for e in range(dim):
                t = t[p[e]]
            logger.info(f'Point {el} added to element {t.i}')
            t.addP(el)
            t.addF(f[I])
            t.addW(self.w[I])
        logger.info('-' * 45)

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
        self.__elemAdd(self.dim, self.elements)

    def AppendLocalMatrix(self, el):
        logger.info('MakeLocalMatrix')
        dim = self.dim
        nums = range(pow(4, dim))
        psi = self.__psi
        local_f_number = 2**dim

        pow_dim_2 = [2**(dim-l) for l in range(dim)]

        L = self.__localMatrixes

        inds = self.indexs
        for I in nums:
            for J in nums:
                value = 1
                for i in range(self.dim):
                    logger.debug(f'L[{i}][{inds[I][i]}][{inds[J][i]}] = {L[i][inds[I][i]][inds[J][i]]}')
                    value *= el.b * L[i][inds[I][i]][inds[J][i]]

                for i, p in enumerate(el.p):
                    value += el.w[i] * psi(el, p, I) * psi(el, p, J)
                
                i = el.nodes[I // local_f_numbers]*local_f_numbers + I % local_f_numbers
                j = el.nodes[J // local_f_numbers]* local_f_numbers + J % local_f_numbers
                logger.debug(f'i={i}\tj={j}')
                self.A[i][j] += value

            for _i, p in enumerate(el.p):
                self.F[i] += el.w[_i] * psi(el, p, I) * el.f[_i]

    def Solve(self):
        logger.info('Solve')
        self.answer = np.linalg.solve(self.A, self.F)
        return self.answer

    def Paint(self):
        logger.info('Paint')
        x = []
        y = []
        z = []
        K = 20
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
            elem_steps = [list(accumulate(np.ones(K, 1) * (self.h[0] / K))), list(accumulate(np.ones(K, 1) * (self.h[1] / K)))]
            for i,el_y in enumerate(self.elements):
                _y = [_el + el_y.mn for _el in elem_steps[1]]
                y.extend(_y)
                for j,el_x in enumerate(el_y):
                    _x = [_el + el_x.mn for _el in elem_steps[0]]
                    if i == 0:
                        x.extend(_x)
                    z.extend([list(accumulate([self.answer[v+i*lfnn] * psi(el, y, v) for v in rle]))[-1] for y in _x])
