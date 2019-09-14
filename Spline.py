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

    def __init__(self, file):
        logger.info('Init')

        self.__f = [Spline.__f1, Spline.__f2, Spline.__f3, Spline.__f4]

        points = np.loadtxt(file)
        logger.info(f'Read {points}')

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

        for el,I in enumerate(points):
            p = [(el[i] - mn[i]) // K for i in range(dim)]
            i = [pow(K,k)*p[k] for k in range(dim)]
            i = accumulate(i)
            logger.info(f'Point {el} added to element {i}')
            self.elements[i].p.append(I)


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
        for I in nums:
            for J in nums:
                value = 0
                for p in el.p:
                    value += w[p] * psi(el, p, I) * psi(el, p, J)
                
    def Solve(self):
        logger.info('Solve')

    def Experiment(self):
        logger.info('Experiment')
