PAL library version 1.5.1.mdw.1

I have made significant changes to the PAL library. These changes are
not authorized by the original authors, and so do not constitute a new
official version of the library. They are based on PAL version 1.5.1

As is the case with the PAL library, these changes are released under
the LGPL version 2.1. (http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)

The following list of changes is not exhaustive. It covers the most
significant changes, and some minor ones, but other minor changes may
have escaped my notice.

======== Substitution Models, Likelihood, MLE Tree Search ==========

* New substitution model substmodel.MosaicSubstitutionModel, which
  allows different substitution models to operate on different parts
  of the tree. The MosaicSubstitutionModel holds a number of
  submodels, and the input tree gets leaves labeled via Annotation to
  say which submodel apply at the leaves.

* Tree edges and nodes can now have additional parameters beyond edge
  length (to support MosaicSubstitutionModel.) The effects of these on
  likelihood calculations depend on the substitution model. These are
  implemented via interfaces substmodel.ModelEdgeParameters and
  substmodel.ModelNodeParameters. There are also 'update' hooks for
  update methods to be called on these as the tree is being created
  (e.g. in
  treesearch.PossiblyRootedMLSearcher.InternalNode.rebuildPattern)

* New interface substmodel.NonTimeReversibleSubstitutionModel, used by
  MosaicSubstitutionModel and new class GeneralNonEquilibriumModel
  which allows any model to have non-equilibrium base frequencies at
  the root.

* Changed treesearch.UnrootedMLSearcher to
  treesearch.PossiblyRootedMLSearcher - can now handle both rooted or
  unrooted tree searches, with non-time-reversible models.(These
  changes are quite extensive.) Substantial comments added explaining
  how this class works internally.

* algorithmics.SearchEngine.run has 'debug' PrintWriter argument to
  which debugging messages are written (set to null for no debugging.)

* New very general rates-across-sites model substmodel.ArbitraryRates

* substmodel.SingleClassSubstitutionModel etc. modified so that the
  equilibrium base frequencies can be considered either as fixed, or
  as model parameters (i.e. when we optimize the model, can base
  frequencies be optimized?)

======== Tree Utilities ==========

* tree.AllRottedTreeIterator and tree.AllUnrootedTreeIterator:
  iterators through all possible tree topologies.

* new class tree.YuleTree for generating random trees by the Yule process.

* new method tree.TreeUtils.stringToTree converts a string containing
  a tree in Newick format into a Tree.

* Bugfixes and upgrades to TreeUtils.reroot. (Turned out to be
  obsolete: TreeManipulator does all this, but nobody had marked
  TreeUtils.reroot as deprecated.)

======== Maths & Maths Utilities ==============

* math.IntMatrix: fast matrix routines when you know the values are
  integer, and won't overflow. (Used by distance.LogDetDistance.)

* math.EigenSystem: Finds eigenvalues/vectors of a matrix, supports
  matrix exponentials and logarithms. Replaces math.MatrixExponential,
  which is buggy (does not correctly handle complex eigenvalues.) Uses
  JAMA linear algebra library.

* PAL now uses JAMA linear algebra library. Uses of the old PAL Matrix class
  should be refactored to use JAMA instead, but this has not been done
  in most cases.

* math.StochasticVector class (mostly for internal use): a
  parameterization of a vector constrained to have nonnegative entries
  and sum to one.

* Pseudorandom numbers generation (math.MersenneTwisterFast): many
  methods changed to allow PRNG seeds to be set. (Otherwise we cannot
  guarantee consistent results between program runs, which among other
  things greatly complicates debugging.) You can cause any use of an
  unseeded PRNG to throw an error with the
  math.MersenneTwisterFast.forbidUnseeded() method.

* class math.SetOfSmallIntegers: implements sets as bitstrings,
  allowing very fast set operations, but a very restricted universal
  set (integers 0 to 63.) Could have its range extended with relative
  ease.


======== Alignments and Distances ============

* Can now perform multiple sequence alignment, via external programs
  (currently only 'muscle'.) See alignment.Align.

* distance.LogDetDistance: calculate logdet distances.


========= Misc changes ===================

* new method misc.Parameterized.ParameterizedUser.setAllParameters

* Minor bug fixes in algorithmics.UndoableAction, math.GammaFunction, 
