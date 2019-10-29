import logging
import numpy as np
import operator
from itertools import accumulate, combinations_with_replacement, permutations

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Net')

class NetFabric:

    inputFile = ''
    kElem = []

    kMx = 1
    kMn = 1

    n = []

    def __init__(self, f, K, kMx=1, kMn=1):
        logger.info('Init')
        self.kElem = np.array(K)
        self.inputFile = f
        self.kMn = kMn
        self.kMx = kMx

    def Generate(self, defW=1, customW={}):
        logger.info('Generate')

        points = np.loadtxt(self.inputFile)
        logger.info(f'Read {points}')

        f = [p[-1] for p in points]
        points = [p[:-1] for p in points]

        mx = []
        dim = len(points[0])
        
        a = [max(points, key=lambda x: x[el])[el] for el in range(dim)]
        mx = np.array(a)
        a = [min(points, key=lambda x: x[el])[el] for el in range(dim)]
        mn = np.array(a)
        K = self.kElem
        mx *= self.kMx
        mn *= self.kMn
        h = (mx - mn) * (1.0 / K)

        w = np.ones(len(points))*defW
        for v in customW:
            for i in customW[v]:
                w[i] = v
        
        N = list(accumulate(np.ones(dim-1) + np.array(K[:-1]), operator.mul))
        N.insert(0, 1)
        self.n = np.array(N)

        kElem = K
        kNode = [k+1 for k in K]
        nElem = list(accumulate(K, operator.mul))[-1]
        nNodes = list(accumulate([el+1 for el in K], operator.mul))[-1]

        self.elements = self.__elemInit(dim-1, [])
        logger.info(f'{nElem} elements created')

        return Net(points,f,)

class Net:

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

    def __init__(self, **kwargs):
        logger.info('Init')
        self.func = func
        self.points = points

        self.elems = kwargs['elems']

        self.elems = kwargs['elems']
        self.elems = kwargs['elems']
        self.elems = kwargs['elems']
        self.elems = kwargs['elems']