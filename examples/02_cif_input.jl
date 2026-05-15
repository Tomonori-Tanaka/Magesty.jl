# Building an SCEBasis from a CIF file (AtomsIO).
#
# Run from the repository root, after installing AtomsIO into the active
# environment (AtomsIO is intentionally NOT a Magesty dependency):
#     julia --project=. -e 'using Pkg; Pkg.add("AtomsIO")'
#     julia --project=. examples/02_cif_input.jl
#
# This example writes a small FePt L1_0 CIF to a temp file, loads it via
# AtomsIO, and feeds the resulting AtomsBase system into the SCEBasis
# constructor. The `interaction` kwarg mirrors the `[interaction]` table
# in the TOML schema.

using Magesty
using AtomsIO

const FEPT_CIF = """
data_FePt_L10
_cell_length_a       2.722
_cell_length_b       2.722
_cell_length_c       3.717
_cell_angle_alpha    90.0
_cell_angle_beta     90.0
_cell_angle_gamma    90.0
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe1 Fe 0.0 0.0 0.0
Pt1 Pt 0.5 0.5 0.5
"""

const INTERACTION = (
    body1 = (lmax = Dict(:Fe => 2, :Pt => 2),),
    body2 = (lsum = 2,
             cutoff = Dict((:Fe, :Fe) => -1.0,
                           (:Fe, :Pt) => -1.0,
                           (:Pt, :Pt) => -1.0)),
)

function main()
    cif_path = tempname() * ".cif"
    write(cif_path, FEPT_CIF)
    try
        system = load_system(cif_path)
        basis = SCEBasis(
            system;
            interaction = INTERACTION,
            tolerance_sym = 1e-5,
            isotropy = false,
            verbosity = false,
        )
        println("Loaded ", length(system), " atoms from CIF.")
        println("Number of SALCs: ", length(basis.salcbasis.salc_list))
        println("Symmetry (international symbol): ", basis.symmetry.international_symbol)
    finally
        isfile(cif_path) && rm(cif_path)
    end
    return nothing
end

main()
