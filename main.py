import logging
import numpy as np
from SplineLine import Spline
# from Element import Element
# from itertools import accumulate
# from math import floor

logging.basicConfig(filename='log.txt', level=logging.INFO)
logger = logging.getLogger('Main')

def GeneratePoints():
    logger.info('GeneratePoints')
    x = np.random.randint(-30, size=(4, 4), high=50)
    np.savetxt('output.txt', x)

def main():
    logger.info('Start')
    # GeneratePoints()

    s = Spline('input.txt')
    s.MakeMatrix()

if __name__ == '__main__':
    main()

