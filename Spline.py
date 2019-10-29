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
from Net import Net

# with open(f'log-{strftime("%H-%M-%S", gmtime())}.txt','w') as f:
#     pass
# logging.basicConfig(filename=f'log-{strftime("%H-%M-%S", gmtime())}.txt', level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Spline')

numberLocalMatrix = 0

class Spline():

    dim = 1

    net = ''

    indexs = []

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
            s *= self.__f[_i](x[k], el.mn[k], self.net.h[k])
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

    def __elemAdd(self, d, mas):
        if d == 0:
            self.AppendLocalMatrix(mas)
            return
        for i in mas:
            self.__elemAdd(d-1, i)

    def __init__(self, net):
        logger.info('Init')

        self.__f = [Spline.__f1, Spline.__f2, Spline.__f3, Spline.__f4]

        self.net = net
        dim = net.dim
        self.dim = dim

        for _h in net.h:
            self.__localMatrixes.append(Spline.__localMatrix(_h))

        l = list(accumulate([k*2 for k in net.kNode], operator.mul))[-1]
        l = range(l)
        self.A = [np.zeros(len(l)) for i in l] 
        self._A = np.copy(self.A)
        self.F = np.zeros(len(l))

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
        self.__elemAdd(self.dim, self.net.elems)

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
