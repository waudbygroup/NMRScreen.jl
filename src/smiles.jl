function smilestoimage(smiles)
    mol = smilestomol(smiles)
    io = IOBuffer()
    drawpng(io, mol, 250, 200) 
    rotr90(load(io))
end
