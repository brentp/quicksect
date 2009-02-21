Quicksect
=========

Description
-----------

Quicksect is a fast python / cython implementation of interval search based on the pure python version in 
`bx-python <http://bx-python.trac.bx.psu.edu/>`__ 

It has a reasonable test-suite in tests/. And the doctests in this
file can be run with

UPDATE: James Taylor has incoporated this code back into bx-python (with his improvements
        and compatibility changes)
        get the code from `bitbucket <http://bitbucket.org/james_taylor/bx-python/>`__



$ nosetests --with-doctest --doctest-extension=.rst README.rst

License is MIT.


Use
---
    >>> from quicksect import IntervalNode, Feature

Feature
+++++++

features are just things with a start, stop, and (optional) strand.

make some fake features...
    >>> starts  = [x for x in range(1, 100, 10)]
    >>> stops   = [s + 4 for s in starts]
    >>> strands = [ i % 2 == 0 and 1 or -1 for i in range(len(starts))]

the feature class:
    >>> feats = [Feature(start, stop, strand) for (start, stop, strand) in zip(starts, stops, strands)]

and potentially a name:
    >>> Feature(1001, 1003, strand=-1, name="fred")
    Feature(1001, 1003, strand=-1, name="fred")



Quicksect
+++++++++

    >>> inter = IntervalNode(feats[0])
    >>> for feat in feats[1:]:
    ...     inter = inter.insert(feat)
    >>> inter                                  #doctest: +ELLIPSIS
    IntervalNode(31, 35)



**search**

    >>> inter.intersect(25, 31)
    [Feature(31, 35, strand=-1), Feature(21, 25, strand=1)]

    >>> inter.intersect(23, 33)
    [Feature(31, 35, strand=-1), Feature(21, 25, strand=1)]



    >>> inter.left(Feature(26, 26), n=1)
    [Feature(21, 25, strand=1)]

    >>> inter.right(Feature(26, 26), n=1)
    [Feature(31, 35, strand=-1)]



