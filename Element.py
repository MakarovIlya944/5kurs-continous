import logging

logger = logging.getLogger('Element')

class Element():

    n = 0

    i = 0
    nodes = []

    f = []
    p = []
    w = []
    b = 0
    mn = 0

    def __init__(self, _mn, _inds):
        self.i = Element.n
        Element.n += 1
        self.mn = _mn
        self.nodes = [int(el) for el in _inds]
        self.p = []
        self.f = []
        self.w = []

    def addP(self, point):
        self.p.append(point)

    def addF(self, point):
        self.f.append(point)

    def addW(self, w):
        self.w.append(w)
