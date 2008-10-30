#!/usr/bin/python
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

#cdef class Interval:
#    cdef public int start, stop
#    def __init__(self, int start, int stop):
#        self.start  = start
#        self.stop   = stop

ctypedef char * char_star

cdef class Feature:
    """\
    Basic feature, with required integer start and stop properties.
    Also accpets optional strand as +1 or -1 (used for up/downstream queries),
    a name, and any arbitrary data is sent in on the info keyword argument

    >>> from quicksect import Feature

    >>> f1 = Feature(23, 36)
    >>> f2 = Feature(34, 48, strand=-1, name="fred", info={'chr':12, 'anno':'transposon'})
    >>> f2
    Feature(34, 48, strand=-1, name="fred", {'anno': 'transposon', 'chr': 12})

    """
    cdef public int start, stop, strand
    cdef public object info
    cdef public object name
    cdef public int chr

    def __init__(self, int start, int stop, int strand=0, int chr=0, object name="", object info=None):
        assert start <= stop, "start must be less than stop"
        self.start  = start
        self.stop   = stop
        self.chr    = chr
        self.strand = strand
        self.name   = name
        self.info   = info

    def __repr__(self):
        fstr = "Feature(%d, %d" % (self.start, self.stop)
        if self.strand != 0:
            fstr += ", strand=%d" % self.strand
        if strlen(self.name) != 0:
            fstr += ', name="' + str(self.name) + '"'
        if not self.info is None:
            fstr += ", " + str(self.info)
        fstr += ")"
        return fstr

cpdef int distance(Feature f1, Feature f2):
    """\
    Distance between 2 features. The integer result is always positive or zero.
    If the features overlap or touch, it is zero.
    >>> from quicksect import Feature, distance
    >>> distance(Feature(1, 2), Feature(12, 13))
    10
    >>> distance(Feature(1, 2), Feature(2, 3))
    0
    >>> distance(Feature(1, 100), Feature(20, 30))
    0

    """
    if f1.stop < f2.start: return f2.start - f1.stop
    if f2.stop < f1.start: return f1.start - f2.stop
    return 0



class IntervalTree:

    def __init__(self):
        self.chrs = {}

    def insert(self, Feature interval):
        chr = interval.chr

        if interval.chr in self.chrs:
            self.chrs[chr] = self.chrs[chr].insert(interval)
        else:
            self.chrs[chr] = IntervalNode(interval)

    def find(self, interval):
        chr = interval.chr
        if chr in self.chrs:
            return self.chrs[chr].intersect(interval.start, interval.stop)

    def traverse( self, func ):
        for item in self.chrs.itervalues():
            item._traverse( func )

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
    >>> from quicksect import IntervalNode, Feature
    >>> tree = IntervalNode(Feature(0, 10, -1))

    Add intervals, the only requirement is that the interval have integer
    start and end attributes. Optional arguments are strand, name, and info.

    >>> Feature(1, 22, strand=-1, name="fred", info={'chr':12, 'anno': 'anything'})
    Feature(1, 22, strand=-1, name="fred", {'anno': 'anything', 'chr': 12})


    >>> tree = tree.insert(Feature(3, 7, 1))
    >>> tree = tree.insert(Feature(3, 40, -1))
    >>> tree = tree.insert(Feature(13, 50, 1))

    Queries
    -------

    find
    ++++

    >>> tree.find(2, 5)
    [Feature(3, 7, strand=1), Feature(3, 40, strand=-1), Feature(0, 10, strand=-1)]
    >>> tree.find(11, 100)
    [Feature(13, 50, strand=1), Feature(3, 40, strand=-1)]
    >>> tree.find(100, 200)
    []

    left/right
    ++++++++++
    the left method finds features that are strictly to the left of
    the query feature. overlapping features are not considered:

    >>> tree.left(Feature(0, 1))
    []
    >>> tree.left(Feature(11, 12))
    [Feature(0, 10, strand=-1)]


    up/downstream
    +++++++++++++
    up/downstream method behave exactly like left/right, except that
    the direction is determined by the strand of the query feature. 
    If the strand is 1, then upstream is left, downstream is right.

    If the strand is -1, then upstream is right, downstream is left.
    >>> tree.upstream(Feature(11, 12, strand=1))
    [Feature(0, 10, strand=-1)]
    >>> tree.upstream(Feature(11, 12, strand=-1))
    [Feature(13, 50, strand=1)]

    all of these method take an argument 'n' for the number of results desired.
    >>> tree.upstream(Feature(1, 2, strand=-1), n=3)
    [Feature(3, 40, strand=-1), Feature(3, 7, strand=1), Feature(13, 50, strand=1)]
    
    nearest neighbors
    +++++++++++++++++
    #>>> tree.nearest_neighbors(Feature(1, 2))
    #[Feature(0, 10, strand=-1)]

    #>>> tree.nearest_neighbors(Feature(1, 2), n=2)
    #[Feature(0, 10, strand=-1), Feature(3, 7, strand=1)]

    """
    cdef float priority
    cdef Feature interval 
    cdef public int start, stop
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
        return "IntervalNode(%i, %i)" % (self.start, self.stop)

    def __cinit__(IntervalNode self, Feature interval):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority   = ceil(nlog * log(-1.0/(1.0 * rand()/RAND_MAX - 1)))
        self.start      = interval.start
        self.stop       = interval.stop
        self.interval   = interval
        self.maxstop    = interval.stop
        self.minstart   = interval.start
        self.minstop    = interval.stop
        self.cleft       = EmptyNode
        self.cright      = EmptyNode
        self.croot       = EmptyNode

    def insert(self, interval):
        return self._insert(interval)

    cdef IntervalNode _insert(IntervalNode self, Feature interval):
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
            self.maxstop = imax3(self.stop, self.cright.maxstop, self.cleft.maxstop)
            self.minstop = imin3(self.stop, self.cright.minstop, self.cleft.minstop)
            self.minstart = imin3(self.start, self.cright.minstart, self.cleft.minstart)
        elif self.cright is not EmptyNode:
            self.maxstop = imax2(self.stop, self.cright.maxstop)
            self.minstop = imin2(self.stop, self.cright.minstop)
            self.minstart = imin2(self.start, self.cright.minstart)
        elif self.cleft is not EmptyNode:
            self.maxstop = imax2(self.stop, self.cleft.maxstop)
            self.minstop = imin2(self.stop, self.cleft.minstop)
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
        #if start <= self.stop and stop >= self.start: results.append(self.interval)
        if (not self.stop < start) and (not self.start > stop): results.append(self.interval)
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

        if -1 < position - self.stop < max_dist:
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

    def neighbors(self, Feature f, int n=1, int max_dist=2500):
        cdef list neighbors = []

        cdef IntervalNode right = self.cright
        while right.cleft is not EmptyNode:
            right = right.cleft

        cdef IntervalNode left = self.cleft
        while left.cright is not EmptyNode:
            left = left.cright
        return [left, right]

    cpdef left(self, Feature f, int n=1, int max_dist=2500):
        """find n features with a start > than f.stop
        f: a Feature object
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use start - 1 becuase .left() assumes strictly left-of
        self._seek_left(f.start - 1, results, n, max_dist)
        if len(results) == n: return results
        r = results
        r.sort(key=operator.attrgetter('stop'), reverse=True)
        return r[:n]

    cpdef right(self, Feature f, int n=1, int max_dist=2500):
        """find n features with a stop < than f.start
        f: a Feature object
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use stop + 1 becuase .right() assumes strictly right-of
        self._seek_right(f.stop + 1, results, n, max_dist)
        if len(results) == n: return results
        r = results
        r.sort(key=operator.attrgetter('start'))
        return r[:n]

    def upstream(self, Feature f, int n=1, int max_dist=2500):
        """find n upstream features where upstream is determined by
        the strand of the query Feature f
        Overlapping features are not considered.

        f: a Feature object
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        if f.strand == -1:
            return self.right(f, n, max_dist)
        return self.left(f, n, max_dist)

    def downstream(self, Feature f, int n=1, int max_dist=2500):
        """find n downstream features where downstream is determined by
        the strand of the query Feature f
        Overlapping features are not considered.

        f: a Feature object
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        if f.strand == -1:
            return self.left(f, n, max_dist)
        return self.right(f, n, max_dist)



    def traverse(self, func):
        self._traverse(func)

    cdef void _traverse(IntervalNode self, object func):
        if self.cleft is not EmptyNode: self.cleft._traverse(func)
        func(self)
        if self.cright is not EmptyNode: self.cright._traverse(func)


cdef IntervalNode EmptyNode = IntervalNode(Feature(0, 0))
