import logging
import numpy as np
from Element import Element
from itertools import accumulate

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
            _i = i % 4
            i /= 4
            # Поправить
            s *= self.__f[_i](x[k], 1, self.h[k])

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

    def __init__(self, file):
        logger.info('Init')

        self.__f = [Spline.__f1, Spline.__f2, Spline.__f3, Spline.__f4]

        points = np.loadtxt(file)
        logger.info(f'Read {points}')

        f = [p[-1] for p in points]
        points = [p[:len(p) - 2] for p in points]

        mx = []
        dim = len(points[0])
        self.dim = dim
        
        a = [max(points, key=lambda x: x[el])[el] for el in range(dim)]
        mx = np.array(a)
        a = [min(points, key=lambda x: x[el])[el] for el in range(dim)]
        mn = np.array(a)
        K = self.K
        self.h = (mx - mn) * (1.0 / K)

        self.elements = [Element() for el in range(pow(K, dim))]
        logger.info(f'{K}^{dim} elements created')

        l = range(4*pow(K, dim))
        self.A = [np.zeros(len(l)) for i in l] 
        self.F = np.zeros(len(l))

        logger.info('-' * 15)
        for I,el in enumerate(points):
            p = [(el[i] - mn[i]) // K for i in range(dim)]
            i = [pow(K,k)*p[k] for k in range(dim)]
            i = list(accumulate(i))[-1]
            logger.info(f'Point {el} added to element {i}')
            self.elements[i].p.append(I)
        logger.info('-' * 45)

        with open(f'{dim}d.txt','r') as f:
            lines = f.readlines()
            self.indexs = [[int(c) for c in str(l)] for l in lines]

    def Calculate(self):
        logger.info('Calculate')
        self.MakeMatrix()
        self.Solve()

    def MakeMatrix(self):
        logger.info('MakeMatrix')
        for el in self.elements:
            self.AppendLocalMatrix(el)
        

    def AppendLocalMatrix(self, el):
        logger.info('MakeLocalMatrix')
        dim = self.dim
        nums = range(pow(4, dim))
        w = self.w
        psi = self.__psi

        L = []
        for h in self.h:
            L.append(Spline.__local_matrix(h))

        for I in nums:
            for J in nums:
                value = 1
                for i in range(self.dim):
                    value *= el.b * L[I*len(nums) + J][i]

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

    def Experiment(self):
        logger.info('Experiment')
