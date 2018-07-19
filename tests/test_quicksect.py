import sys, os
import unittest
try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
except:
    sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

from quicksect import Interval, IntervalNode, IntervalTree, distance

import quicksect
try:
    from cPickle import dumps, loads
except ImportError:
    from pickle import dumps, loads
import operator
import random

class NeighborTestCase(unittest.TestCase):

    def setUp(self):
        iv = IntervalNode(Interval(50, 59))
        for i in range(0, 110, 10):
            if i == 50: continue
            f = Interval(i, i + 9)
            iv = iv.insert(f)
        self.intervals = iv

    def test_left(self):
        iv = self.intervals 
        self.assertEqual(str(iv.left(Interval(60, 70), n=2)), str([Interval(50, 59), Interval(40, 49)]))

        for i in range(10, 100, 10):
            f = Interval(i, i)
            r = iv.left(f, max_dist=10, n=1)
            self.assertEqual(r[0].end,  i - 1)

    def test_toomany(self):
        iv = self.intervals
        self.assertEqual(len(iv.left(Interval(60, 70), n=200)) , 6)


    def test_right(self):
        iv = self.intervals
        self.assertEqual(str(iv.left(Interval(60, 70), n=2)), str([Interval(50, 59), Interval(40, 49)]))

        def get_right_start(b10):
            r = iv.right(Interval(b10, b10 + 1), n=1)
            assert len(r) == 1
            return r[0].start
        
        for i in range(10, 100, 10):
            self.assertEqual(get_right_start(i), i + 10)

        for i in range(0, 100, 10):
            f = Interval(i - 1, i - 1)
            r = iv.right(f, max_dist=10, n=1)
            self.assertEqual(r[0].start, i)


    def test_n(self):
        iv = self.intervals
        for i in range(0, 90, 10):
            f = Interval(i + 1, i + 1)
            r = iv.right(f, max_dist=20, n=2)
            self.assertEqual(r[0].start, i + 10)
            self.assertEqual(r[1].start, i + 20)


class RelativeTestCase(unittest.TestCase):
    def setUp(self):
        intervals = []
        for i in range(11, 20000, 15):
            for zz in range(random.randint(2, 5)):
                m = random.randint(1, 10)
                p = random.randint(1, 10)
                intervals.append(Interval(i - m, i + p))
        iv = IntervalNode(intervals[0])
        for f in intervals[1:]:
            iv = iv.insert(f)

        self.intervals = intervals
        self.tree = iv

    def test_left(self):
        max_dist = 200
        n = 15
        iv = self.tree
        for i in range(11, 20000, 25):
            for zz in range(random.randint(2, 5)):
                s1 = random.randint(i + 1, i + 20 )
                f = Interval(s1, s1)

                bf = brute_force_find_left(self.intervals, f, max_dist, n)
                tf = iv.left(f, max_dist=max_dist, n=n)
                if len(tf) == 0:
                    assert len(bf) == 0, bf
                    continue


                mdist = max(distance(f, t) for t in tf)
                self.assertTrue(set(bf).issuperset(tf))
                diff = set(bf).difference(tf)
                self.assertTrue(len(diff) == 0, (diff))


    def test_right(self):
        max_dist = 200
        n = 15
        iv = self.tree
        for i in range(11, 20000, 25):
            for zz in range(random.randint(1, 6)):
                s1 = random.randint(i + 1, i + 20 )
                f = Interval(s1, s1)

                bf = brute_force_find_right(self.intervals, f, max_dist, n)
                tf = iv.right(f, max_dist=max_dist, n=n)
                if len(tf) == 0:
                    assert len(bf) == 0, bf
                    continue


                mdist = max(distance(f, t) for t in tf)
                self.assertTrue(set(bf).issuperset(tf))
                diff = set(bf).difference(tf)
                self.assertTrue(len(diff) == 0, (diff))


class LotsaTestCase(unittest.TestCase):
    """ put lotsa data in the tree and make sure it works"""
    def setUp(self):
        iv = IntervalNode(Interval(1, 2))
        self.max = 1000000
        for i in range(0, self.max, 10):
            f = Interval(i, i)
            iv = iv.insert(f)

        for i in range(6000):
            iv = iv.insert(Interval(0, 1))
        self.intervals = iv



    def test_count(self):
        iv = self.intervals
        
        r = iv.right(Interval(1, 1), n=33)
        self.assertEqual(len(r), 33)

        l = iv.left(Interval(1, 1), n=33)
        self.assertEqual(len(l), 1)


    def test_max_dist(self):
        iv = self.intervals
        r = iv.right(Interval(1, 1), max_dist=0, n=10)
        self.assertEqual(len(r), 0)

        for n, d in enumerate(range(10, 1000, 10)):
            r = iv.right(Interval(1, 1), max_dist=d, n=10000)
            self.assertEqual(len(r), n + 1) 

    def test_find(self):
        iv = self.intervals

        intervals = []
        iv.traverse(intervals.append)

        for t in range(250):
            start = random.randint(0, self.max - 10000)
            stop  = start + random.randint(100, 10000)

            results = iv.find(start, stop)
            for feat in results:
                self.assertTrue(
                        (feat.end >= start and feat.end <= stop) 
                            or 
                        (feat.start <= stop and feat.start >= start)
                        )
            bf = brute_force_find(intervals, start, stop)
            assert len(results) == len(bf)

def brute_force_find(intervals, start, stop):
    return [i for i in intervals if i.end >= start and i.start <= stop]

def brute_force_find_left(intervals, f, max_dist, n):
    r = [x for x in brute_force_find(intervals, 0, f.start)\
               if x.end < f.start and distance(x, f) <= max_dist]
    r.sort(key=operator.attrgetter('end'), reverse=True)
    if len(r) <= n: return r
    i = n
    while distance(r[i], f) == distance(r[i - 1], f):
        i += 1
    return r[:i]



def brute_force_find_right(intervals, f, max_dist, n):
    r = [x for x in brute_force_find(intervals, f.end, 99999999999)\
               if x.start > f.start and distance(x, f) <= max_dist]
    r.sort(key=operator.attrgetter('start'))
    if len(r) <= n: return r
    i = n
    while distance(r[i], f) == distance(r[i - 1], f):
        i += 1
    return r[:i]


class PickleTestCase(unittest.TestCase):
    """ test pickling."""
    def setUp(self):
        pass

    def test_feature_pickle(self):
        f = Interval(22, 38, data={'a':22})
        g = loads(dumps(f))
        self.assertEqual(f.start, g.start)
        self.assertEqual(g.data['a'], 22)

    def test_tree_pickle(self):
        a = IntervalTree()
        for ichr in range(5):
            for i in range(10, 100, 6):
                f = Interval(i -4, i + 4)
                a.insert(f)
        
        a.dump('a.pkl')

        b = IntervalTree()
        b.load('a.pkl')
        for ichr in range(5):
            for i in range(10, 100, 6):
                f = Interval(i -4, i + 4)
                af = sorted(a.find(f), key=operator.attrgetter('start'))
                bf = sorted(b.find(f), key=operator.attrgetter('start'))

                assert len(bf) > 0
                self.assertEqual(len(af), len(bf))
                self.assertEqual(af[0].start, bf[0].start)
                self.assertEqual(af[-1].start, bf[-1].start)

class EmptyTreeTestCase(unittest.TestCase):
    """ test search on an empty tree."""

    def setUp(self):
        self.tree = IntervalTree()

    def test_search(self):
        self.tree.search(46, 47)

    def test_find(self):
        self.tree.find(Interval(46, 47))

    def test_left(self):
        self.tree.left(Interval(46, 47))

    def test_right(self):
        self.tree.right(Interval(46, 47))


class TestIssue9(unittest.TestCase):
    def setUp(self):
        self.tree4 = IntervalTree()
        self.tree4.insert(Interval(22, 33, data='example1'))
        self.tree4.insert(Interval(22, 33, data='example2'))

    def test_right(self):
        self.assertEqual(0, len(self.tree4.right(Interval(44,55))))
        self.assertEqual(2, len(self.tree4.right(Interval(11,12))))

    def test_left(self):
        self.assertEqual(2, len(self.tree4.left(Interval(44,55))))
        self.assertEqual(0, len(self.tree4.left(Interval(11,12))))

def main():
    unittest.main()



if __name__ == "__main__":
    unittest.main()
