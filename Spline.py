import logging
import numpy as np

logger = logging.getLogger('Spline')

class Spline():

    def __init__(self):
        logger.info('Init')
        x = self.ReadPoints('input.txt')
        self.MakeElements(x)

    def MakeElements(self, points):
        logger.info('MakeElements')

        mxx = max(points, key=lambda x: x[0])[0]
        mxy = max(points, key=lambda x: x[1])[1]
        mxz = max(points, key=lambda x: x[2])[2]

        mnx = min(points, key=lambda x: x[0])[0]
        mny = min(points, key=lambda x: x[1])[1]
        mnz = min(points, key=lambda x: x[2])[2]

        mx = np.array([mxx,mxy,mxz])
        mn = np.array([mnx,mny,mnz])
        dx = mx - mn
        k = 10

        


    def ReadPoints(self, file):
        logger.info('ReadPoints')
        x = np.loadtxt(file)
        logger.info(f'Read {x}')
        return x
