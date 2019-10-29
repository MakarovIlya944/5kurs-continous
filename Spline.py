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

# with open(f'log-{strftime("%H-%M-%S", gmtime())}.txt','w') as f:
#     pass
# logging.basicConfig(filename=f'log-{strftime("%H-%M-%S", gmtime())}.txt', level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Spline')

numberLocalMatrix = 0

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
            for i in range(self.K[self.dim-d-1]):
                mn = list(np.append(mas, i)*self.h)
                mn = np.array(mn)
                res.append(Element(mn + self.mn, self.__elem(np.append(mas, i))))
            return res
        return [self.__elemInit(d-1, np.append(mas, i)) for i in range(self.K[self.dim-d-1])]

    def __elemAdd(self, d, mas):
        if d == 0:
            self.AppendLocalMatrix(mas)
            return
        for i in mas:
            self.__elemAdd(d-1, i)

    def __init__(self, file, _K, w=None):
        logger.info('Init')

        self.__f = [Spline.__f1, Spline.__f2, Spline.__f3, Spline.__f4]
        self.K = np.array(_K)

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
        k_part = 1.0010
        mx *= k_part
        self.h = (mx - mn) * (1.0 / K)
        h = self.h

        self.w = np.ones(len(points))
        # self.w[5]=0
        
        for _h in self.h:
            self.__localMatrixes.append(Spline.__localMatrix(_h))

        N = list(accumulate(np.ones(dim-1) + np.array(K[:-1]), operator.mul))
        N.insert(0, 1)
        self.n = np.array(N)

        self.kElem = K
        self.kNode = [k+1 for k in K]
        kn = self.kNode
        self.nElem = list(accumulate(K, operator.mul))[-1]
        self.nNodes = list(accumulate([el+1 for el in K], operator.mul))[-1]

        self.elements = self.__elemInit(dim-1, [])
        logger.info(f'{self.nElem}^{dim} elements created')

        l = list(accumulate([k*2 for k in kn], operator.mul))[-1]
        l = range(l)
        self.A = [np.zeros(len(l)) for i in l] 
        self._A = np.copy(self.A)
        self.F = np.zeros(len(l))

        logger.debug('-' * 45)
        for I,el in enumerate(points):
            p = np.floor((el-mn)/h)
            p = list([int(i) if K[ind] != int(i) else K[ind] - 1 for ind,i in enumerate(p)])
            t = self.elements
            for e in range(dim):
                t = t[p[e]]
            logger.debug(f'Point {el} added to element {t.i}')
            t.addP(el)
            t.addF(f[I])
            t.addW(self.w[I])
        logger.debug('-' * 45)

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
                
                i = el.nodes[I // local_f_number] * local_f_number + I % local_f_number
                j = el.nodes[J // local_f_number] * local_f_number + J % local_f_number
                logger.debug(f'i={i}\tj={j}')
                
                if logger.level == logging.DEBUG:
                    self._A[i][j] = numberLocalMatrix + 1
                    # self._A[i][j] = value
                self.A[i][j] += value

            for _i, p in enumerate(el.p):
                self.F[i] += el.w[_i] * psi(el, p, I) * el.f[_i]
        if logger.level == logging.DEBUG:
            np.savetxt(f'step_{numberLocalMatrix}.txt',self._A,fmt='%.0f')
            numberLocalMatrix += 1

    def Solve(self):
        logger.info('Solve')
        self.answer = np.linalg.lstsq(self.A, self.F, rcond=None)
        return self.answer
