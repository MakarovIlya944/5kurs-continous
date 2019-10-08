import argparse
import logging
import numpy as np
from Spline import Spline
import matplotlib.pyplot as plt

logging.basicConfig(filename='log.txt', level=logging.WARNING)
logger = logging.getLogger('Main')

def func(x):
    return x*x

def GeneratePoints(noise=True, dim=1, f=func):
    logger.info('GeneratePoints')
    k = 10
    if dim == 1:
        x = np.random.randint(-k, size=(10, 1), high=k)
        x = [np.array([i,func(i) + (x[i] if noise else 0) * (1 - abs(i - 5)/5)]) for i in range(0,10)]
    elif dim == 2:
        x = []
        for i in range(-5,5):
            for j in range(-5,5):
                x.append(np.array([i, j, func(i) + (np.random.randint(-k, size=(1,1), high=k) if noise else 0)]))
    np.savetxt('input.txt', x)

def main():
    logger.info('Start')
    GeneratePoints(False, 2)

    s = Spline('input.txt', 2, 5)
    s.MakeMatrix()
    np.savetxt('before_solveA.txt',s.A, fmt='%1.2f')
    np.savetxt('before_solveF.txt',s.F, fmt='%1.2f')
    # plt.matshow(s.A)
    ans = s.Solve()
    np.savetxt('answer.txt',ans, fmt='%1.2f')


    s.Paint()

if __name__ == '__main__':
    main()
