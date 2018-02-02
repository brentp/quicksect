#!/usr/bin/python2.5
"""
Intersects ... faster.  Suports GenomicInterval datatype and multiple
chromosomes.
"""
import operator

#@cdef extern from "math.h":
cdef extern from "stdlib.h":
    int ceil(float f)
    float log(float f)
    int RAND_MAX
    int rand()
    int strlen(char *)
    int iabs(int)

cdef class Interval:
    cdef public int start, end
    cdef public object data
    def __init__(self, int start, int end, data=None):
        self.start = start
        self.end = end
        self.data = data
    def __repr__(self):
        if self.data is not None:
            return "Interval(%d, %d, data=%s)" % (self.start, self.end, self.data)
        else:
            return "Interval(%d, %d)" % (self.start, self.end)
    def __reduce__(self):
        args = self.__getstate__()
        return type(self), (args.pop('start'), args.pop('end')), args

    def __getstate__(self):
        return {'start': self.start, 'end': self.end, 'data': self.data}

    def __setstate__(self, kwargs):
        self.data = kwargs['data']

cpdef int distance(Interval f1, Interval f2):
    """\
    Distance between 2 features. The integer result is always positive or zero.
    If the features overlap or touch, it is zero.
    >>> from quicksect import Interval, distance
    >>> distance(Interval(1, 2), Interval(12, 13))
    10
    >>> distance(Interval(1, 2), Interval(2, 3))
    0
    >>> distance(Interval(1, 100), Interval(20, 30))
    0

    """
    if f1.end < f2.start: return f2.start - f1.end
    if f2.end < f1.start: return f1.start - f2.end
    return 0


cdef class IntervalTree:
    cdef IntervalNode root

    def __init__(self):
        self.root = None

    cpdef insert(self, Interval interval):

        if self.root is None:
            self.root = IntervalNode(interval)
        else:
            self.root = self.root.insert(interval)

    def add(self, int start, int end, other=None):
        return self.insert(Interval(start, end, other))

    def find(self, interval):
        if self.root is None:
            return []
        else:
            return self.root.intersect(interval.start, interval.end)

    def search(self, int start, int end):
        if self.root is None:
            return []
        else:
            return self.root.intersect(start, end)

    def left(self, Interval f, int n=1, int max_dist=25000):
        if self.root is None:
            return []
        else:
            return self.root.left(f, n, max_dist)

    def right(self, Interval f, int n=1, int max_dist=25000):
        if self.root is None:
            return []
        else:
            return self.root.right(f, n, max_dist)

    def dump(self, fn):
        try:
            import cPickle
        except ImportError:
            import pickle as cPickle
        l = []
        a = l.append
        self.root.traverse(a)
        fh = open(fn, "wb")
        for f in l:
            cPickle.dump(f, fh)

    def load(self, fn):
        try:
            import cPickle
        except ImportError:
            import pickle as cPickle
        fh = open(fn, "rb")
        while True:
            try:
                feature = cPickle.load(fh)
                self.insert(feature)
            except EOFError:
                break


cdef inline int imax2(int a, int b):
    if b > a: return b
    return a

cdef inline int imax3(int a, int b, int c):
    if b > a:
        if c > b:
            return c
        return b
    if a > c:
        return a
    return c

cdef inline int imin3(int a, int b, int c):
    if b < a:
        if c < b:
            return c
        return b
    if a < c:
        return a
    return c

cdef inline int imin2(int a, int b):
    if b < a: return b
    return a

cdef float nlog = -1.0 / log(0.5)

cdef class IntervalNode:
    """\
    Data structure for performing intersect and neighbor queries on a 
    set of intervals. Algorithm uses a segment/interval tree to perform
    efficient queries. 

    Usage
    =====
    >>> from quicksect import IntervalNode, Interval
    >>> tree = IntervalNode(Interval(0, 10))

    Add intervals, the only requirement is that the interval have integer
    start and end attributes. Optional arguments are strand, name, and info.

    >>> Interval(1, 22, info={'chr':12, 'anno': 'anything'})


    >>> tree = tree.insert(Interval(3, 7, 1))
    >>> tree = tree.insert(Interval(3, 40, -1))
    >>> tree = tree.insert(Interval(13, 50, 1))

    Queries
    -------

    find
    ++++

    >>> tree.find(2, 5)
    [Interval(3, 7), Interval(3, 40), Interval(0, 10)]
    >>> tree.find(11, 100)
    [Interval(13, 50), Interval(3, 40)]
    >>> tree.find(100, 200)
    []

    left/right
    ++++++++++
    the left method finds features that are strictly to the left of
    the query feature. overlapping features are not considered:

    >>> tree.left(Interval(0, 1))
    []
    >>> tree.left(Interval(11, 12))
    [Interval(0, 10)]

    """
    cdef float priority
    cdef public Interval interval
    cdef public int start, end
    cdef int minstop, maxstop, minstart
    cdef IntervalNode cleft, cright, croot

    property left_node:
        def __get__(self):
            return self.cleft if self.cleft is not EmptyNode else None
    property right_node:
        def __get__(self):
            return self.cright if self.cright is not EmptyNode else None
    property root_node:
        def __get__(self):
            return self.croot if self.croot is not EmptyNode else None
    


    def __repr__(self):
        return "IntervalNode(%i, %i)" % (self.start, self.end)

    def __cinit__(IntervalNode self, Interval interval):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority   = ceil(nlog * log(-1.0/(1.0 * rand()/RAND_MAX - 1)))
        self.start      = interval.start
        self.end       = interval.end
        self.interval   = interval
        self.maxstop    = interval.end
        self.minstart   = interval.start
        self.minstop    = interval.end
        self.cleft       = EmptyNode
        self.cright      = EmptyNode
        self.croot       = EmptyNode

    def insert(self, interval):
        return self._insert(interval)

    cdef IntervalNode _insert(IntervalNode self, Interval interval):
        cdef IntervalNode croot = self
        if interval.start > self.start:

            # insert to cright tree
            if self.cright is not EmptyNode:
                self.cright = self.cright._insert(interval )
            else:
                self.cright = IntervalNode(interval)
            # rebalance tree
            if self.priority < self.cright.priority:
                croot = self.rotate_left()

        else:
            # insert to cleft tree
            if self.cleft is not EmptyNode:
                self.cleft = self.cleft._insert(interval)
            else:
                self.cleft = IntervalNode(interval)
            # rebalance tree
            if self.priority < self.cleft.priority:
                croot = self.rotate_right()
    
        croot.set_stops()
        self.cleft.croot  = croot
        self.cright.croot = croot
        return croot

    cdef IntervalNode rotate_right(IntervalNode self):
        cdef IntervalNode croot = self.cleft
        self.cleft  = self.cleft.cright
        croot.cright = self
        self.set_stops()
        return croot

    cdef IntervalNode rotate_left(IntervalNode self):
        cdef IntervalNode croot = self.cright
        self.cright = self.cright.cleft
        croot.cleft  = self
        self.set_stops()
        return croot

    cdef inline void set_stops(IntervalNode self):
        if self.cright is not EmptyNode and self.cleft is not EmptyNode: 
            self.maxstop = imax3(self.end, self.cright.maxstop, self.cleft.maxstop)
            self.minstop = imin3(self.end, self.cright.minstop, self.cleft.minstop)
            self.minstart = imin3(self.start, self.cright.minstart, self.cleft.minstart)
        elif self.cright is not EmptyNode:
            self.maxstop = imax2(self.end, self.cright.maxstop)
            self.minstop = imin2(self.end, self.cright.minstop)
            self.minstart = imin2(self.start, self.cright.minstart)
        elif self.cleft is not EmptyNode:
            self.maxstop = imax2(self.end, self.cleft.maxstop)
            self.minstop = imin2(self.end, self.cleft.minstop)
            self.minstart = imin2(self.start, self.cleft.minstart)
        

    def intersect(self, int start, int stop):
        """
        given a start and a stop, return a list of features
        falling within that range
        """
        cdef list results = []
        self._intersect(start, stop, results)
        return results

    find = intersect
        
    cdef void _intersect(IntervalNode self, int start, int stop, list results):
        # to have starts, stops be non-inclusive, replace <= with <  and >= with >
        #if start <= self.end and stop >= self.start: results.append(self.interval)
        if (not self.end < start) and (not self.start > stop): results.append(self.interval)
        #if self.cleft is not EmptyNode and start <= self.cleft.maxstop:
        if self.cleft is not EmptyNode and not self.cleft.maxstop < start:
            self.cleft._intersect(start, stop, results)
        #if self.cright is not EmptyNode and stop >= self.start:
        if self.cright is not EmptyNode and not self.start > stop:
            self.cright._intersect(start, stop, results)
    

    cdef void _seek_left(IntervalNode self, int position, list results, int n, int max_dist):
        # we know we can bail in these 2 cases.
        if self.maxstop + max_dist < position: return
        if self.minstart > position: return

        #import sys;sys.stderr.write( " ".join(map(str, ["SEEK_LEFT:", self, self.cleft, self.maxstop, self.minstart,  position])))

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cright is not EmptyNode:
                self.cright._seek_left(position, results, n, max_dist)

        if -1 < position - self.end < max_dist:
            results.append(self.interval)

        # TODO: can these conditionals be more stringent?
        if self.cleft is not EmptyNode:
                self.cleft._seek_left(position, results, n, max_dist)


    
    cdef void _seek_right(IntervalNode self, int position, list results, int n, int max_dist):
        # we know we can bail in these 2 cases.
        if self.maxstop < position: return
        if self.minstart - max_dist > position: return

        #print "SEEK_RIGHT:",self, self.cleft, self.maxstop, self.minstart, position

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cleft is not EmptyNode: 
                self.cleft._seek_right(position, results, n, max_dist)

        if -1 < self.start - position < max_dist:
            results.append(self.interval)

        if self.cright is not EmptyNode:
                self.cright._seek_right(position, results, n, max_dist)

    def neighbors(self, Interval f, int n=1, int max_dist=25000):
        cdef list neighbors = []

        cdef IntervalNode right = self.cright
        while right.cleft is not EmptyNode:
            right = right.cleft

        cdef IntervalNode left = self.cleft
        while left.cright is not EmptyNode:
            left = left.cright
        return [left, right]

    cpdef left(self, Interval f, int n=1, int max_dist=25000):
        """find n features with a start > than f.end
        f: a Interval object
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use start - 1 becuase .left() assumes strictly left-of
        self._seek_left(f.start - 1, results, n, max_dist)
        if len(results) <= n: return results
        r = results
        r.sort(key=operator.attrgetter('end'), reverse=True)
        if distance(f, r[n]) != distance(f, r[n-1]):
            return r[:n]
        while distance(r[n], f) == distance(r[n - 1], f):
            n += 1
        return r[:n]


    cpdef right(self, Interval f, int n=1, int max_dist=25000):
        """find n features with a stop < than f.start
        f: a Interval object
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use stop + 1 becuase .right() assumes strictly right-of
        self._seek_right(f.end + 1, results, n, max_dist)
        if len(results) <= n: return results
        r = results
        r.sort(key=operator.attrgetter('start'))
        if distance(f, r[n]) != distance(f, r[n-1]):
            return r[:n]
        while distance(r[n], f) == distance(r[n - 1], f):
            n += 1
        return r[:n]

    def __iter__(self):
            
        if self.cleft is not EmptyNode:
            yield self.cleft

        yield self.interval

        if self.cright is not EmptyNode:
            yield self.cright

    def traverse(self, func):
        self._traverse(func)


    cdef void _traverse(IntervalNode self, object func):
        if self.cleft is not EmptyNode: self.cleft._traverse(func)
        func(self.interval)
        if self.cright is not EmptyNode: self.cright._traverse(func)


cdef IntervalNode EmptyNode = IntervalNode(Interval(0, 0))
