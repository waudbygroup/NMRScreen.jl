function smilestoimage(smiles, size=(250,200))
    mol = smilestomol(smiles)
    io = IOBuffer()
    drawpng(io, mol, size...) 
    rotr90(load(io))
end
