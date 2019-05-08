Fortran snippets
=========================================

Some algorithms implemented in Fortran.
This is not a library but a collection of modules to be re-used in other codes.
The algorithms must be:
 1. General enough to be re-used
 2. No external dependencies (to depend on other modules inside snippets is allowed)
 3. Each snippet should include the module and some examples or benchmarks

hashtable
--------------------------
The hashtable can be used to implement dictionaries or sets.
The example implmentation is used to map integers to integers.
It should be possible to use for other data structures by choosing a different
hash function and minor changes in the code.

octree
--------------------------
Octree is a data structure used among others to find the nearest neighbour of a point
in a list of other points in a three dimensional space. It works by creating a big box
that contains all the space and then recursively spliting it in octets and assigning the
points to each of them.
It is possible then to query for the nerest neighbour very efficientely.
This structure can be used for example to map 3D points of one list to another.

quadratic
--------------------------
Routines to fit a quandratic function in 3D, interpolate it and integrate
using an hybrid tetrahedron method.
