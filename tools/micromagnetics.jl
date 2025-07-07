"""
Caculation tool for the derivation of the micromagnetics model parameters
"""
using ArgParse
using LinearAlgebra
using EzXML
using JLD2


function calc_micromagnetics(input_xml::String, input_system::String, cutoff::Float64)::Nothing
    # Read the XML file
    doc = readxml(input_xml)

    # Read the system jld2 file
    system = load(input_system)

    sce_basis_set = findfirst("//SCEBasisSet", doc)
    if isnothing(sce_basis_set)
        throw(ArgumentError("<SCEBasisSet> node not found in the XML file."))
    end

    for salc in EzXML.findall("SALC", sce_basis_set)
        index = parse(Int, salc["index"])
    end

end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_xml"
        help = "Input xml file"
        required = true

        "--input_system"
        help = "Input system jld2 file"
        required = true

        "--cutoff", "-c"
        help = "Cutoff radius"
        default = -1.0
    end

    args = parse_args(s)
end

main()