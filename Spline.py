import logging
import numpy as np
from Element import Element
from itertools import accumulate
from math import floor

logger = logging.getLogger('Spline')

class Spline():

    w = []
    a = []
    elements = []
    points = []
    K = 10
    dim = 1
    h = []
    local_matrix = np.array([ 
        np.array([12, 6, -12, 6]),
        np.array([6, 4, -6, 2]),
        np.array([-12, -6, 12, -6]),
        np.array([6, 2, -6, 4])])
    indexs = []

    A = []
    F = []

    __f = []

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

    def __psi(self, el, p, i):
        s = 1
        x = self.points[p]
        for k in range(self.dim):
            _i = self.indexs[i][k]
            s *= self.__f[_i](x[k], el.mn[k], self.h[k])

    def __local_matrix(h):
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

    def __elemInit(d):
        if d == 0:
            return 

    def __init__(self, file):
        logger.info('Init')

        self.__f = [Spline.__f1, Spline.__f2, Spline.__f3, Spline.__f4]

        points = np.loadtxt(file)
        logger.info(f'Read {points}')

        f = [p[-1] for p in points]
        a = len(points[0]) - 1
        points = [p[:a] for p in points]

        mx = []
        dim = len(points[0])
        self.dim = dim
        
        a = [max(points, key=lambda x: x[el])[el] for el in range(dim)]
        mx = np.array(a)
        a = [min(points, key=lambda x: x[el])[el] for el in range(dim)]
        mn = np.array(a)
        K = self.K
        self.h = (mx - mn) * (1.0 / K)
        h = self.h
        k_part = 1.01
        mx *= k_part

        # self.elements = [Element() for el in range(pow(K, dim))]
        rrange = range(K)
        result = []
        for d in dim:
            result = [Element() for z in rrange]
            

        for x in rrange:
            tmp = []
            for y in rrange:
                tmp.append([Element() for z in rrange])
            self.elements.append(tmp)

        logger.info(f'{K}^{dim} elements created')

        l = range(4*pow(K, dim))
        self.A = [np.zeros(len(l)) for i in l] 
        self.F = np.zeros(len(l))

        logger.info('-' * 15)
        for I,el in enumerate(points):
            p = [floor((el[i] - mn[i]) / h[i]) for i in range(dim)]
            p = [i if K != i else K - 1 for i in p]
            t = self.elements
            for e in range(dim):
                t = t[p[e]]
            logger.info(f'Point {el} added to element {t.i}')
            t.p.append(el)
            t.f.append(f[I])
        logger.info('-' * 45)

        with open(f'{dim}d.txt','r') as f:
            lines = f.readlines()
            # Not implemented dim>10
            self.indexs = [[int(c) for c in str(int(l))] for l in lines]

    def Calculate(self):
        logger.info('Calculate')
        self.MakeMatrix()
        self.Solve()

    def MakeMatrix(self):
        logger.info('MakeMatrix')
        for x in self.elements:
            for y in x:
                for z in y:
                    self.AppendLocalMatrix(z)
        

    def AppendLocalMatrix(self, el):
        logger.info('MakeLocalMatrix')
        dim = self.dim
        nums = range(pow(4, dim))
        w = self.w
        psi = self.__psi

        L = []
        for h in self.h:
            L.append(Spline.__local_matrix(h))

        inds = self.indexs
        for I in nums:
            for J in nums:
                value = 1
                K = I*len(nums) + J
                for i in range(self.dim):
                    value *= el.b * L[i][inds[K][I]][inds[K][J]]

                for p in el.p:
                    value += w[p] * psi(el, p, I) * psi(el, p, J)

                i = el.p[I // 4] * 4 + I % 4
                j = el.p[J // 4] * 4 + J % 4
                self.A[i][j] += value

            for p in el.p:
                self.F[i] += w[p] * psi(el, p, I)

    def Solve(self):
        logger.info('Solve')
        return np.linalg.solve(self.A, self.F)