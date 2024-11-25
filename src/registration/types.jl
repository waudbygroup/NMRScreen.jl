struct RegistrationCocktail
    id
    name
    fragment_ids
    smiles
    refspec
    boundspec
    peak_ids
    library_shifts
    ref_shifts
    bound_shifts
    good::Observable{Vector{Bool}}
end
