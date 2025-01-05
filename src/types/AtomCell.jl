module AtomCells

import Base: isless, ==, hash, show

export AtomCell
"""
	AtomCell(atom::Int, cell::Int)
A structure that stores atom index and imaginary (virtual) cell index
"""
struct AtomCell
	atom::Int
	cell::Int

	function AtomCell(atom::Integer, cell::Integer)
		if (cell < 1) || (cell > 27)
			throw(ArgumentError("cell: $cell. the cell should be 1 <= cell <= 27."))
		end
		new(atom, cell)
	end
end

isless(ac1::AtomCell, ac2::AtomCell) = (ac1.atom, ac1.cell) < (ac2.atom, ac2.cell)
isless(ac1::AtomCell, tuple::NTuple{2, Integer}) = (ac1.atom, ac1.cell) < tuple
isless(tuple::NTuple{2, Integer}, ac2::AtomCell) = tuple < (ac2.atom, ac2.cell)
==(ac1::AtomCell, ac2::AtomCell) = (ac1.atom, ac1.cell) == (ac2.atom, ac2.cell)
==(ac1::AtomCell, tuple::NTuple{2, Integer}) = (ac1.atom, ac1.cell) == tuple
==(tuple::NTuple{2, Integer}, ac2::AtomCell) = tuple == (ac2.atom, ac2.cell)
hash(ac::AtomCell, h::UInt) = hash((ac.atom, ac.cell), h)
show(io::IO, ac::AtomCell) = print(io, "(atom: $(ac.atom), cell: $(ac.cell))")

end
