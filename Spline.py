from Element import Element
from itertools import accumulate, combinations_with_replacement, permutations
from math import floor
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
import operator
from mpl_toolkits.mplot3d import Axes3D
from time import gmtime, strftime
import numpy as np
import logging

with open(f'log-{strftime("%H-%M-%S", gmtime())}.txt','w') as f:
    pass
logging.basicConfig(filename=f'log-{strftime("%H-%M-%S", gmtime())}.txt', level=logging.DEBUG)
logger = logging.getLogger('Main')

