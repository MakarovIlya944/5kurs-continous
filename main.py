import argparse
import logging
import numpy as np
from Spline import Spline
import matplotlib.pyplot as plt

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
        np.savetxt('input.txt', points)
        return data
            
def main():
    logger.info('Start')

    f = PointsFabric(1, lambda x: x*x,lambda x: 0,[[-1,1,0.2]])
    potins = f.generate()

    s = Spline('input.txt', [2], 5)
    s.MakeMatrix()
    np.savetxt('before_solveA.txt',s.A, fmt='%1.2e')
    np.savetxt('before_solveF.txt',s.F, fmt='%1.2f')
    # plt.matshow(s.A)
    ans = s.Solve()
    np.savetxt('answer.txt',ans, fmt='%1.2f')
    s.Paint(potins)

if __name__ == '__main__':
    main()
