import logging
import numpy as np
from Spline import Spline
# from Element import Element
# from itertools import accumulate, combinations_with_replacement, permutations
# from math import floor
# import matplotlib.pyplot as plt
# from matplotlib import cm
# import operator
# from mpl_toolkits.mplot3d import Axes3D

logging.basicConfig(filename='log.txt', level=logging.DEBUG)
logger = logging.getLogger('Main')

def GeneratePoints(random=True):
    logger.info('GeneratePoints')
    if random:
        x = np.random.randint(0, size=(20, 3), high=10)
    else:
        k = 10
        x = np.random.randint(-k, size=(10, 1), high=k)
        x = [np.array([i,i*i + x[i] * (1 - abs(i - 5)/5)]) for i in range(0,10)]
    np.savetxt('input.txt', x)

def main():
    logger.info('Start')
    GeneratePoints()

    s = Spline('input.txt')
    s.MakeMatrix()
    np.savetxt('before_solveA.txt',s.A, fmt='%1.2e')
    np.savetxt('before_solveF.txt',s.F, fmt='%1.2e')
    ans = s.Solve()
    np.savetxt('answer.txt',ans)
    s.Paint()

if __name__ == '__main__':
    main()

