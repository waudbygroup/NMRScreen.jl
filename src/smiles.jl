function smilestoimage(smiles)
    mol = smilestomol(smiles)
    io = IOBuffer()
    drawpng(io, mol, 300, 200) 
    rotr90(load(io))
end
