import sys, os
import unittest
try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
except:
    sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

from quicksect import Feature
from quicksect import IntervalNode

import quicksect

class NeighborTestCase(unittest.TestCase):

    def setUp(self):
        iv = IntervalNode(Feature(50, 59))
        for i in range(0, 110, 10):
            if i == 50: continue
            f = Feature(i, i + 9)
            iv = iv.insert(f)
        self.intervals = iv

    def test_left(self):
        iv = self.intervals 
        self.assertEqual(str(iv.left(Feature(60, 70), n=2)), str([Feature(50, 59), Feature(40, 49)]))

        for i in range(10, 100, 10):
            f = Feature(i, i)
            r = iv.left(f, max_dist=10, n=1)
            self.assertEqual(r[0].stop,  i - 1)

    def test_toomany(self):
        iv = self.intervals 
        self.assertEqual(len(iv.left(Feature(60, 70), n=200)) , 6)


    def test_right(self):
        iv = self.intervals 
        self.assertEqual(str(iv.left(Feature(60, 70), n=2)), str([Feature(50, 59), Feature(40, 49)]))

        def get_right_start(b10):
            r = iv.right(Feature(b10, b10 + 1), n=1)
            assert len(r) == 1
            return r[0].start
        
        for i in range(10, 100, 10):
            self.assertEqual(get_right_start(i), i + 10)

        for i in range(0, 100, 10):
            f = Feature(i - 1, i - 1)
            r = iv.right(f, max_dist=10, n=1)
            self.assertEqual(r[0].start, i)

    def test_upstream(self):
        iv = self.intervals 
        upstreams = iv.upstream(Feature(59, 60), n=200)
        for u in upstreams:
            self.assertTrue(u.stop < 59)

        upstreams = iv.upstream(Feature(60, 70, strand=-1), n=200)
        for u in upstreams:
            self.assertTrue(u.start > 70)


        upstreams = iv.upstream(Feature(58, 58, strand=-1), n=200)
        for u in upstreams:
            self.assertTrue(u.start > 59)


    def test_downstream(self):
        iv = self.intervals 
        downstreams = iv.downstream(Feature(59, 60), n=200)
        for d in downstreams:
            self.assertTrue(d.start > 60)

        downstreams = iv.downstream(Feature(59, 60, strand=-1), n=200)
        for d in downstreams:
            self.assertTrue(d.start < 59)
        

    def test_n(self):
        iv = self.intervals
        for i in range(0, 90, 10):
            f = Feature(i + 1, i + 1)
            r = iv.right(f, max_dist=20, n=2)
            self.assertEqual(r[0].start, i + 10)
            self.assertEqual(r[1].start, i + 20)


class LotsaTestCase(unittest.TestCase):
    """ put lotsa data in the tree and make sure it works"""
    def setUp(self):
        iv = IntervalNode(Feature(1, 2))
        self.max = 1000000
        for i in range(0, self.max, 10):
            f = Feature(i, i)
            iv = iv.insert(f)

        for i in range(600):
            iv = iv.insert(Feature(0, 1))
        self.intervals = iv



    def test_count(self):
        iv = self.intervals
        
        r = iv.right(Feature(1, 1), n=33)
        self.assertEqual(len(r), 33)

        l = iv.left(Feature(1, 1), n=33)
        self.assertEqual(len(l), 1)

        u = iv.upstream(Feature(1, 1, strand=-1), n=9999)
        self.assertEqual(len(u), 250)

        # now increase max_dist
        u = iv.upstream(Feature(1, 1, strand=-1), n=9999, max_dist=99999)
        self.assertEqual(len(u), 9999)


    def test_max_dist(self):        
        iv = self.intervals
        r = iv.right(Feature(1, 1), max_dist=0, n=10)
        self.assertEqual(len(r), 0)

        for n, d in enumerate(range(10, 1000, 10)):
            r = iv.right(Feature(1, 1), max_dist=d, n=10000)
            self.assertEqual(len(r), n + 1) 

    def test_find(self):
        iv = self.intervals
        import random
        for t in range(25):
            start = random.randint(0, self.max - 10000)
            stop  = start + random.randint(100, 10000)

            results = iv.find(start, stop)
            for feat in results:
                self.assertTrue(
                        (feat.stop >= start and feat.stop <= stop) 
                            or 
                        (feat.start <= stop and feat.start >= start)
                        )

    #def test_iter(self):
    #    tree = self.intervals
    #    for item in tree:
    #        print item

from cPickle import dumps, loads
class PickleTestCase(unittest.TestCase):
    """ test pickling."""
    def setUp(self):
        pass

    def test_feature_pickle(self):
        f = Feature(22, 38, strand=-1, chr='ddd',name="fred", info={'a':22})
        g = loads(dumps(f))
        self.assertEqual(f.start, g.start)
        self.assertEqual(g.info['a'], 22)


if __name__ == "__main__":
    unittest.main()
