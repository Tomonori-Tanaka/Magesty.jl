"""
Caculation tool for the derivation of the micromagnetics model parameters
"""
using ArgParse
using LinearAlgebra
using EzXML
using JLD2

include("../src/Structures.jl")
using .Structures


function calc_micromagnetics(input_xml::String, cutoff::Float64)::Nothing
    # Read the XML file
    doc = readxml(input_xml)

    structure = Structure(input_xml)
    println(structure)


    sce_basis_set = findfirst("//SCEBasisSet", doc)
    if isnothing(sce_basis_set)
        throw(ArgumentError("<SCEBasisSet> node not found in the XML file."))
    end
    # Get the number of SALCs
    num_salc = parse(Int, sce_basis_set["num_salc"])

    for i_salc in 1:num_salc
        salc_node = findfirst("//SALC[index='$i_salc']", sce_basis_set)
    end

end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_xml"
        help = "Input xml file"
        required = true

        "--cutoff", "-c"
        help = "Cutoff radius"
        default = -1.0
    end

    args = parse_args(s)
end

main()