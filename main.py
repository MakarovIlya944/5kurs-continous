import logging
import numpy as np
# from Spline import Spline
from Element import Element
from itertools import accumulate, combinations_with_replacement, permutations
from math import floor
import matplotlib.pyplot as plt
from matplotlib import cm
import operator
from mpl_toolkits.mplot3d import Axes3D

logging.basicConfig(filename='log.txt', level=logging.WARNING)
logger = logging.getLogger('Main')

def GeneratePoints(random=True, dim=1):
    logger.info('GeneratePoints')
    if random:
        x = np.random.randint(0, size=(20, 3), high=10)
    else:
        if dim == 1:
            k = 10
            x = np.random.randint(-k, size=(10, 1), high=k)
            x = [np.array([i,i*i + x[i] * (1 - abs(i - 5)/5)]) for i in range(0,10)]
        elif dim == 2:
            x = []
            for i in range(-5,5):
                for j in range(-5,5):
                    x.append(np.array([i, j, i*i]))
    np.savetxt('input.txt', x)

def main():
    logger.info('Start')
    GeneratePoints(False, 2)

    s = Spline('input.txt')
    s.MakeMatrix()
    np.savetxt('before_solveA.txt',s.A, fmt='%1.2e')
    np.savetxt('before_solveF.txt',s.F, fmt='%1.2f')
    ans = s.Solve()
    np.savetxt('answer.txt',ans, fmt='%1.2f')
    

    s.Paint()

if __name__ == '__main__':
    main()

