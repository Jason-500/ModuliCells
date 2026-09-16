from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.topology.cell_complex import GenericCellComplex


class ModuliCellComplex(GenericCellComplex):
    r"""
    Finite cell complex equipped with finite cell-symmetry data.

    The underlying cells are stored dimension by dimension.  The class is
    intentionally agnostic about the combinatorial model used to describe a
    cell: a cell may be represented by a DFV tuple, a ribbon graph, or any
    future object exposing an ``automorphism_group()`` method.

    ``global_isotropy_order`` records a finite subgroup acting trivially on
    every effective combinatorial cell.  For the genus-two Weierstrass model
    this is the hyperelliptic involution, of order 2.
    """

    def __init__(self, cells, boundary_cells=None, global_isotropy_order=1):
        if not cells:
            self._cells = {}
            self._dimension = -1
        else:
            self._dimension = max(cells)
            self._cells = {
                d: set(cells.get(d, ()))
                for d in range(self._dimension + 1)
            }

        if boundary_cells is None:
            self._boundary_cells = {
                d: set() for d in range(self._dimension + 1)
            }
        else:
            self._boundary_cells = {
                d: set(boundary_cells.get(d, ()))
                for d in range(self._dimension + 1)
            }

        global_isotropy_order = ZZ(global_isotropy_order)
        if global_isotropy_order <= 0:
            raise ValueError("global_isotropy_order must be positive")
        self._global_isotropy_order = global_isotropy_order

    def _repr_(self):
        if self._dimension < 0:
            return "Empty moduli cell complex"
        return (
            "Moduli cell complex of dimension {} with {} cells".format(
                self._dimension,
                sum(len(cells) for cells in self._cells.values()),
            )
        )

    def dimension(self):
        return self._dimension

    def cells(self, subcomplex=None):
        r"""
        Return the cells, indexed by dimension.

        If ``subcomplex`` is supplied, remove the cells occurring in the
        subcomplex, following the convention of ``GenericCellComplex``.
        """
        if subcomplex is None:
            return {d: set(cells) for d, cells in self._cells.items()}

        subcells = subcomplex.cells()
        return {
            d: set(cells).difference(subcells.get(d, ()))
            for d, cells in self._cells.items()
        }

    def boundary_cells(self, dimension=None):
        r"""
        Return cells classified as degeneration/bordification cells.

        These cells are metadata attached to the constructed complex; they
        are not included in ``cells()`` and hence do not contribute to the
        ordinary Euler characteristic of the open moduli-space complex.
        """
        if dimension is None:
            return {
                d: set(cells) for d, cells in self._boundary_cells.items()
            }
        return set(self._boundary_cells.get(dimension, ()))

    def global_isotropy_order(self):
        return self._global_isotropy_order

    def effective_stabilizer_order(self, cell):
        r"""Return the order of the effective combinatorial symmetry group."""
        try:
            group = cell.automorphism_group()
        except AttributeError as err:
            raise TypeError(
                "cells must implement automorphism_group() in order to "
                "compute orbifold Euler characteristics"
            ) from err
        return ZZ(group.order())

    def stabilizer_order(self, cell):
        r"""
        Return the full stabilizer order used for orbifold weighting.

        The current model factors the stabilizer as an effective symmetry
        group of the combinatorial cell times a global ineffective isotropy
        group.  This method is deliberately isolated so that more general
        stacky stabilizer models can override it later.
        """
        return self._global_isotropy_order * self.effective_stabilizer_order(cell)

    def effective_orbifold_euler_characteristic(self):
        r"""Return the symmetry-weighted Euler characteristic before global isotropy."""
        chi = QQ.zero()
        for d, cells in self._cells.items():
            for cell in cells:
                chi += QQ((-1) ** d) / self.effective_stabilizer_order(cell)
        return chi

    def orbifold_euler_characteristic(self):
        r"""Return the orbifold Euler characteristic of the stored quotient model."""
        chi = QQ.zero()
        for d, cells in self._cells.items():
            for cell in cells:
                chi += QQ((-1) ** d) / self.stabilizer_order(cell)
        return chi
