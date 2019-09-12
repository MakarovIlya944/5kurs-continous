import logging
import numpy as np
from Spline import Spline

logger = logging.getLogger('Main')

def GeneratePoints():
    logger.info('GeneratePoints')
    x = np.random.randint(-3, size=(4, 4), high=5)
    np.savetxt('output.txt', x)


def main():
    logger.info('Start')
    # GeneratePoints()

    s = Spline()
    s.ReadPoints('input.txt')

if __name__ == '__main__':
    main()

