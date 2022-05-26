function PFATT(spectra_target, initial_conc, X, S, k)
    # Selected spectra
    spectra_target = spectra_target[:]
    spectra = X[:, spectra_target]

    # Estimate final concentration
    initial_conc = initial_conc[:]
    final_conc = 1 - sum(initial_conc, 2)

    # Concentration matrix
    conc_matrix = [initial_conc final_conc]
    conc_matrix = pinv(conc_matrix)

    # Row matrix
    Row_matrix = S.matrixResuts.u * S.matrixResuts.s

    # Diagonal matrix
    Diag_matrix = Diagonal(inv.(S.errorsResults.explvar[1:k]))

    # Transformation matrix
    tmp1 = spectra * conc_matrix
    Transformation_matrix = Diag_matrix * Row_matrix' * tmp1

    # Real factor matrix
    Real_factors_matrix = Row_matrix * Transformation_matrix

    # Real concentration matrix
    Real_conc_matrix = inv(Transformation_matrix) * S.matrixResuts.V'

    # Real spectra reproduction
    Real_spectra_reprod = Real_factors_matrix * Real_conc_matrix

    # Error in reproduction
    rms = norm(S.dr .- Real_spectra_reprod) / norm(S.dr)

    # Wrap up results
    TransformationResults(Row_matrix, Diag_matrix, Transformation_matrix, Real_factors_matrix,
                          Real_conc_matrix, Real_spectra_reprod, rms)
end
