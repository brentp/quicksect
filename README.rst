Quicksect
=========

Description
-----------

Quicksect is a fast python / cython implementation of interval search based on the pure python version in 
`bx-python <http://bx-python.trac.bx.psu.edu/>`__ 
I pulled it out, optimized and converted to cython and James Taylor has incoporated it back into bx-python
with his improvements.

I have brought this project back from the dead because I want a fast, simple, no-dependencies Interval
tree.


$ python setup.py test

License is MIT.

Use
---
    >>> from quicksect import IntervalNode, Interval, IntervalTree

Most common use will be via IntervalTree:

    >>> tree = IntervalTree()
    >>> tree.add(23, 45)
    >>> tree.add(55, 66)
    >>> tree.search(46, 47)
    []
    >>> tree.search(44, 56)
    [Interval(55, 66), Interval(23, 45)]

    >>> tree.insert(Interval(88, 444))
    >>> res = tree.find(Interval(99, 100))
    >>> res
    [Interval(88, 444)]
    >>> res[0].start, res[0].end
    (88, 444)

Thats pretty much everything you need to know about the tree.

Low-Level
+++++++++

In some cases, users may want to utilize the lower-level interface that accesses
the nodes of the tree:

    >>> inter = IntervalNode(Interval(22, 33))
    >>> inter = inter.insert(Interval(44, 55))
    >>> inter.intersect(24, 26)
    [Interval(22, 33)]

    >>> inter.left(Interval(34, 35), n=1)
    [Interval(22, 33)]

    >>> inter.right(Interval(34, 35), n=1)
    [Interval(44, 55)]
