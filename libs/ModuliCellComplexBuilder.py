from libs.CustomErrors import DimensionMismatchError, NotWellDefinedError
from libs.ModuliCellComplex import ModuliCellComplex


class ModuliCellComplexBuilder:
    """Construct a moduli cell complex recursively from top-dimensional DFV cells."""

    def __init__(self, initial_dfv_list):
        if not initial_dfv_list:
            raise ValueError("at least one top-dimensional cell is required")

        self.top_cell_dim = initial_dfv_list[0].dim

        for dfv in initial_dfv_list:
            if dfv.dim != self.top_cell_dim:
                raise DimensionMismatchError(
                    "initial cells must all have the same dimension"
                )
            if not dfv.is_canonical:
                raise NotWellDefinedError(
                    "initial cells must be canonical representatives"
                )
            if not dfv.is_interior:
                raise NotWellDefinedError(
                    "initial cells must lie in the interior"
                )
            if not dfv.is_root_cell:
                raise NotWellDefinedError(
                    "top-dimensional cells must be marked as root cells"
                )

        self.interior_cell_dict = {
            self.top_cell_dim: list(initial_dfv_list)
        }
        self.boundary_cell_dict = {
            self.top_cell_dim: []
        }
        self.boundary_dim = self.top_cell_dim

    @staticmethod
    def _canonicalize(face, canonical_cells):
        """Attach ``face`` to an existing canonical representative if possible."""
        for canonical_cell in canonical_cells:
            is_isomorphic, permutation = (
                canonical_cell.is_orientation_preserving_isomorphic_to(face)
            )
            if is_isomorphic:
                face.is_canonical = False
                face.canonical_image = (canonical_cell, permutation)
                return canonical_cell

        face.set_canonical()
        canonical_cells.append(face)
        return face

    def solve_boundary(self, solve_boundary_cells=False):
        """Construct and canonicalize all facets one dimension lower."""
        if self.boundary_dim <= 0:
            return

        current_cells = self.interior_cell_dict[self.boundary_dim]
        new_interior_cells = []
        new_boundary_cells = []

        for dfv in current_cells:
            interior_faces, boundary_faces = dfv.facets(
                solve_boundary_cells=solve_boundary_cells
            )

            for face in interior_faces:
                self._canonicalize(face, new_interior_cells)

            if solve_boundary_cells:
                for face in boundary_faces:
                    self._canonicalize(face, new_boundary_cells)

        self.boundary_dim -= 1
        self.interior_cell_dict[self.boundary_dim] = new_interior_cells
        self.boundary_cell_dict[self.boundary_dim] = new_boundary_cells

    def solve_boundary_till_dim(self, dim, solve_boundary_cells=False):
        """Recursively construct cells until dimension ``dim`` is reached."""
        if dim < 0 or dim > self.boundary_dim:
            raise ValueError(
                "target dimension must satisfy 0 <= dim <= current dimension"
            )

        while self.boundary_dim != dim:
            self.solve_boundary(
                solve_boundary_cells=solve_boundary_cells
            )

    def cell_complex(self, global_isotropy_order=1):
        """Return the completed result as a ``ModuliCellComplex``."""
        if self.boundary_dim != 0:
            raise NotWellDefinedError(
                "the decomposition is incomplete; call "
                "solve_boundary_till_dim(0) first"
            )

        return ModuliCellComplex(
            self.interior_cell_dict,
            boundary_cells=self.boundary_cell_dict,
            global_isotropy_order=global_isotropy_order,
        )
