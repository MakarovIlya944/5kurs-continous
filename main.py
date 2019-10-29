import argparse
import logging
import numpy as np
from Spline import Spline
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

    def generate(self):
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

    dim = 1
    if dim == 1:
        f = lambda x: x*x
        q = lambda x: (5-abs(x-5)*1)*2
        domains = [[0,10,1]]
        elements = [2]
    elif dim == 2:
        f = lambda x,y: x*x
        q = lambda x,y: 0
        domains = [[0,1,0.1],[0,1,0.1]]
        elements = [1,1]
    elif dim == 3:
        f = lambda x,y,z: (x*x+y*y)*z
        q = lambda x,y,z: ((x*x+y*y)*z)*np.random.rand(1)*0.01 + (100 if abs(x-2.5) < 1 and abs(y-2.5) < 1 else 0)
        domains = [[0,5,0.5],[0,5,0.5],[-1,1,0.2]]
        elements = [2,2,1]

    K = 10

    f = PointsFabric(dim, f, q, domains)
    clear, noise = f.generate()

    s = Spline('input.txt', elements, K)
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
            if len(a) == s.nNodes * (2**dim):
                np.savetxt(f'data/answer_{i}.txt', a, fmt='%1.2f')
                i += 1
        except Exception:
            pass

    p = Painter('data/answer_0.txt', dim, True, s._Spline__psi, K, s.elements, s.mx, s.kElem, s.h, clearPoints=noise)

    testFunc = lambda x,y,z: x+z
    p.Paint(True)

if __name__ == '__main__':
    main()
