import argparse
import logging
import numpy as np
from Spline import Spline
from Net import NetFabric
from Painter import Painter
import matplotlib.pyplot as plt
import math

logging.basicConfig(filename='log.txt', level=logging.WARNING)
logger = logging.getLogger('Main')

class PointsFabric():

    funFile = 'functions.txt'
    # main function
    f = ''
    # noise
    q = ''
    # area of functions
    r = ''
    d = 1
    rnd = False

    def __init__(self, dim, f, q, r=[[0,5,1]], random=False):
        self.f = f
        self.q = q
        self.r = [np.arange(r[i][0],r[i][1],r[i][2]) for i in range(dim)] 
        self.d = dim
        self.rnd = random

    def Generate(self):
        points = []
        data = []
        if self.d == 1:
            for x in self.r[0]:
                points.append(np.array([x, self.f(x) + self.q(x)]))
                data.append(np.array([x, self.f(x)]))
        elif self.d == 2:
            for x in self.r[0]:
                for y in self.r[1]:
                    points.append(np.array([x, y, self.f(x,y) + self.q(x,y)]))
                    data.append(np.array([x, y, self.f(x,y)]))
        elif self.d == 3:
            for x in self.r[0]:
                for y in self.r[1]:
                    for z in self.r[2]:
                        points.append(np.array([x, y, z, self.f(x,y,z) + self.q(x,y,z)]))
                        data.append(np.array([x, y, z, self.f(x,y,z)]))
        np.savetxt('input.txt', points)
        return data, points
  
def main():
    logger.info('Start')

    isShowMatrix = False
    isSaveMatrix = False

    dim = 3
    if dim == 1:
        f = lambda x: x*x
        q = lambda x: (5-abs(x-5)*1)*np.random.rand(1) + (10 if abs(x-5) < 1 else 0) 
        domains = [[0,10,0.5]]
        elements = [3]
    elif dim == 2:
        f = lambda x,y: x*x + y*y
        q = lambda x,y: 100 if abs(x-2) < 1 and abs(y-2) < 1 else 0
        domains = [[-5,5,1],[-5,5,1]]
        elements = [2,1]
    elif dim == 3:
        f = lambda x,y,z: (x*x+y*y)*z
        q = lambda x,y,z: ((x*x+y*y)*z)*np.random.rand(1)*0.01 + (100 if abs(x-2) < 1 and abs(y-2) < 1 else 0)
        domains = [[0,5,1],[0,5,1],[-1,1,0.2]]
        elements = [2,2,1]

    K = 10

    f = PointsFabric(dim, f, q, domains)
    clear, noise = f.Generate()

    n = NetFabric('input.txt', elements, kMn=1.01)
    n = n.Generate(customW={0:[i for i in range(121,130)]})#customB={0:[[0,0]]}, customW={0:[10,11]}

    s = Spline(n)
    s.MakeMatrix()

    if isSaveMatrix:
        np.savetxt('data/before_solveA.txt', s.A, fmt='%1.2e')
        np.savetxt('data/before_solveF.txt', s.F, fmt='%1.2f')
    if isShowMatrix:
        plt.matshow(s.A)
    ans = s.Solve()

    i = 0
    for a in ans:
        try:
            if len(a) == n.nNode * (2**dim):
                np.savetxt(f'data/answer_{i}.txt', a, fmt='%1.2f')
                i += 1
        except Exception:
            pass

    p = Painter('data/answer_0.txt', dim, True, s._Spline__psi, K, n.elems, n.mx, n.kElem, n.h, clearPoints=clear, noisePoints=noise)
    p.Paint(True)

if __name__ == '__main__':
    main()
