import logging
import numpy as np
# from SplineLine import Spline
from Element import Element
from itertools import accumulate
from math import floor
import matplotlib.pyplot as plt

logging.basicConfig(filename='log.txt', level=logging.DEBUG)
logger = logging.getLogger('Main')

def GeneratePoints(random=True):
    logger.info('GeneratePoints')
    if random:
        x = np.random.randint(-30, size=(4, 2), high=50)
    else:
        x = [np.array([i,i*i]) for i in range(0,25)]
    np.savetxt('input.txt', x)

def main():
    logger.info('Start')
    GeneratePoints(False)

    s = Spline('input.txt')
    s.MakeMatrix()
    np.savetxt('before_solveA.txt',s.A)
    np.savetxt('before_solveF.txt',s.F)
    ans = s.Solve()
    np.savetxt('answer.txt',ans)
    s.Paint()

if __name__ == '__main__':
    main()

