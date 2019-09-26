import logging
import numpy as np
from Element import Element
from itertools import accumulate
import matplotlib.pyplot as plt
from math import floor

logger = logging.getLogger('Spline')

class Spline():

    w = []
    a = []
    elements = []
    points = []
    K = 3
    dim = 1
    h = []
    local_matrix = np.array([ 
        np.array([12, 6, -12, 6]),
        np.array([6, 4, -6, 2]),
        np.array([-12, -6, 12, -6]),
        np.array([6, 2, -6, 4])])

    mn = []
    mx = []
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

    def __psi(self, el, x, i):
        s = 1
        for k in range(self.dim):
            _i = self.indexs[i][k]
            s *= self.__f[_i](x[k], el.mn[k], self.h[k])
        return s

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

    def __elemInit(self, d, mas):
        if d == 0:
            return [Element(np.append(mas, i)*self.h + self.mn, np.append(mas, i)) for i in range(self.K)]
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

        # self.elements = [Element() for el in range(pow(K, dim))]
        self.elements = self.__elemInit(dim-1, [])
        logger.info(f'{K}^{dim} elements created')

        l = range((K+1)*pow(2, dim))
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
            # t.p.append(el)
            # t.f.append(f[I])
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

        L = []
        for h in self.h:
            L.append(Spline.__local_matrix(h))

        inds = self.indexs
        for I in nums:
            for J in nums:
                value = 1
                # K = I*len(nums) + J
                for i in range(self.dim):
                    logger.debug(f'L[{i}][{inds[I][i]}][{inds[J][i]}] = {L[i][inds[I][i]][inds[J][i]]}')
                    value *= el.b * L[i][inds[I][i]][inds[J][i]]

                for i, p in enumerate(el.p):
                    value += el.w[i] * psi(el, p, I) * psi(el, p, J)

                i = list(accumulate([((I // local_f_number)%(2**(l+1))+el.indexes[l])*(2**(dim-l)) for l in range(dim)]))[-1] + I % local_f_number
                j = list(accumulate([((J // local_f_number)%(2**(l+1))+el.indexes[l])*(2**(dim-l)) for l in range(dim)]))[-1] + J % local_f_number
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
        # function = []
        # dirivate = []
        y = []
        K = 5
        psi = self.__psi
        if self.dim == 1:
            elem_steps = [self.h[0] * (el/K) for el in range(K)]
            for i,el in enumerate(self.elements):
                _x = [_el + el.mn for _el in elem_steps]
                x.extend(_x)
                # function.extend([self.answer[0+i*4] * psi(el, y, 0+i*4) + self.answer[2+i*4] * psi(el, y, 2+i*4) for y in _x])
                # dirivate.extend([self.answer[1+i*4] * psi(el, y, 1+i*4) + self.answer[3+i*4] * psi(el, y, 3+i*4) for y in _x])
                y.extend([list(accumulate([self.answer[v+i*2] * psi(el, y, v) for v in range(4)]))[-1] for y in _x])
            plt.plot(x,y,'-',self.points, self.f, 'o')
            plt.show()