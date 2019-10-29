import logging
from Element import Element
import numpy as np
import operator
from itertools import accumulate, combinations_with_replacement, permutations

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Net')

class NetFabric:

    inputFile = ''
    kElem = []

    dim = 1

    kMx = 1
    kMn = 1

    b = 0
    n = []
    mn = []
    h = []

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
            for i in range(self.kElem[self.dim-d-1]):
                mn = list(np.append(mas, i)*self.h)
                mn = np.array(mn)
                res.append(Element(mn + self.mn, self.__elem(np.append(mas, i)),b=self.b))
            return res
        return [self.__elemInit(d-1, np.append(mas, i)) for i in range(self.kElem[self.dim-d-1])]

    def __init__(self, f, K, kMx=1, kMn=1):
        logger.info('Builder init')
        self.kElem = np.array(K)
        self.inputFile = f
        self.kMn = kMn
        self.kMx = kMx

    def Generate(self, defW=1, customW={}, defB=0, customB={}):
        logger.info('Generate')

        points = np.loadtxt(self.inputFile)
        logger.info(f'Read {len(points)} points')
        logger.debug(str(points))

        args = {}

        f = [p[-1] for p in points]
        args['f'] = f
        points = [p[:-1] for p in points]
        args['points'] = points

        mx = []
        dim = len(points[0])
        self.dim = dim
        args['dim'] = dim
        
        a = [max(points, key=lambda x: x[el])[el] for el in range(dim)]
        mx = np.array(a)
        a = [min(points, key=lambda x: x[el])[el] for el in range(dim)]
        mn = np.array(a)
        K = self.kElem
        mx *= self.kMx
        mn *= self.kMn
        h = (mx - mn) * (1.0 / K)
        self.mn = mn
        self.h = h
        args['mx'] = mx
        args['mn'] = mn
        args['h'] = h

        w = np.ones(len(points))*defW
        for v in customW:
            logger.debug(f'w={v}, i:{customW[v]}')
            for i in customW[v]:
                w[i] = v
        args['w'] = w

        N = list(accumulate(np.ones(dim-1) + np.array(K[:-1]), operator.mul))
        N.insert(0, 1)
        self.n = np.array(N)

        args['kElem'] = K
        args['kNode'] = [k+1 for k in K]
        nElem = list(accumulate(K, operator.mul))[-1]
        args['nElem'] = nElem
        args['nNodes'] = list(accumulate([el+1 for el in K], operator.mul))[-1]

        self.b = defB
        args['elems'] = self.__elemInit(dim-1, [])
        logger.info(f'{nElem} elements created')

        for I,el in enumerate(points):
            p = np.floor((el-mn)/h)
            p = list([int(i) if K[ind] != int(i) else K[ind] - 1 for ind,i in enumerate(p)])
            t = args['elems']
            for e in range(dim):
                t = t[p[e]]
            logger.debug(f'Point {el} added to element {t.i}')
            t.addP(el)
            t.addF(f[I])
            t.addW(w[I])

        for b in customB:
            logger.debug(f'b={b}, i:{customB[b]}')
            for i in customB[b]:
                t = args['elems']
                for e in range(dim):
                    t = t[i[e]]
                t.b = b

        return Net(kwargs=args)

class Net:

    dim = 1

    elems = []
    h = []

    mn = []
    mx = []

    nNode = 0
    nElem = 0
    kNode = []
    kElem = []

    w = []
    points = []
    func = []

    def __init__(self, kwargs):
        logger.info('Init')
        self.elems = kwargs['elems']
        self.h = kwargs['h']
        self.mn = kwargs['mn']
        self.mx = kwargs['mx']
        self.nNode = kwargs['nNodes']
        self.nElem = kwargs['nElem']
        self.kNode = kwargs['kNode']
        self.kElem = kwargs['kElem']
        self.w = kwargs['w']
        self.points = kwargs['points']
        self.func = kwargs['f']
        self.dim = kwargs['dim']