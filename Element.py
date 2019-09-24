import logging

logger = logging.getLogger('Element')

class Element():

    f = []
    p = []
    b = 1
    i = 0
    n = 0
    mn = 0

    def __init__(self):
        self.i = Element.n
        Element.n += 1