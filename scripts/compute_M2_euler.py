from sage.all import QQ

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from data.DFVdata import D_raw, F_raw, V_raw
from libs.DFV import DFV
from libs.Functions import poly
from libs.ModuliCellComplexBuilder import ModuliCellComplexBuilder


def genus2_top_cells():
    """Return the canonical six-dimensional interior DFV cells for M_2."""
    indices = []

    for i in range(len(D_raw)):
        parameter_polyhedron = poly(D_raw[i], F_raw[i])
        center = parameter_polyhedron.center()

        if (
            parameter_polyhedron.dim() == 6
            and all(x != 0 and x != 1 for x in center)
        ):
            indices.append(i)

    return [
        DFV(
            D_raw[i],
            F_raw[i],
            V_raw[i],
            is_interior=True,
            is_canonical=True,
            is_root_cell=True,
        )
        for i in indices
    ]


def build_genus2_complex(solve_boundary_cells=False):
    """Construct the genus-two cell decomposition from its top cells."""
    top_cells = genus2_top_cells()
    builder = ModuliCellComplexBuilder(top_cells)
    builder.solve_boundary_till_dim(
        0,
        solve_boundary_cells=solve_boundary_cells,
    )

    # The effective DFV action misses the universal hyperelliptic
    # involution, which acts trivially on the Weierstrass quotient data.
    return builder.cell_complex(global_isotropy_order=2)


def main():
    complex_m2 = build_genus2_complex()

    print("number of top cells =", len(complex_m2.n_cells(6)))

    chi_eff = QQ.zero()
    for d in range(complex_m2.dimension() + 1):
        cells = complex_m2.n_cells(d)
        weighted_number = sum(
            (QQ.one() / complex_m2.effective_stabilizer_order(cell)
             for cell in cells),
            QQ.zero(),
        )
        contribution = (-1) ** d * weighted_number
        chi_eff += contribution

        print(
            "dim =", d,
            "number =", len(cells),
            "weighted =", weighted_number,
            "contribution =", contribution,
        )

    print("effective orbifold Euler characteristic =", chi_eff)
    print(
        "orbifold Euler characteristic of M_2 =",
        complex_m2.orbifold_euler_characteristic(),
    )


if __name__ == "__main__":
    main()
