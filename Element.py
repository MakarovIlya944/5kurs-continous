import logging

logger = logging.getLogger('Element')

class Elemets():

    elems = []

    def Init(self, Elements):
        

class Element():

    count = 0
    Id = 0
    PointsId = []

    def __init__(self):
        logger.info('Init')
        self.Id = Element.count
        Element.count += 1
