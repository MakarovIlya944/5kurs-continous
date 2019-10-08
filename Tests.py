from main import *
import unittest

class TestOneDim(unittest.TestCase):
    def FullAnswer(self):
        print()

class TestTwoDim(unittest.TestCase):
    def a(self):
        print()

class TestThreeDim(unittest.TestCase):
    def a(self):
        print()

def test():
    print('Test')

    splineTestSuite = unittest.TestSuite()

    splineTestSuite.addTest(unittest.makeSuite(TestOneDim))
    splineTestSuite.addTest(unittest.makeSuite(TestTwoDim))
    splineTestSuite.addTest(unittest.makeSuite(TestThreeDim))

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(splineTestSuite)
    if len(result.errors) == 0 or len(result.failures) == 0:
        exit(0)
    else:
        exit(1)

if __name__=='__main__':
    test()