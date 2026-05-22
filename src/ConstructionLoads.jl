


struct ConstructionLoadInputs

    P
    W1
    W2
    W3
    L

    E
    I

    aMn       # per-ft capacity (lb·in/ft)
    aVn       # per-ft capacity (lb/ft)

    unit_width  # strip width in ft — used to normalise strip-total demands to per-ft for interaction check

end


struct ConstructionLoadOutputs

    inputs

    M_plus_1
    M_plus_2
    M_plus_3
    M_plus

    M_neg_1
    M_neg_2
    M_neg

    P_ext_1
    P_ext_2
    P_ext_3

    P_int_1
    P_int_2
    P_int_3

    Δ

    interaction

end


function simple_span_M(p)

    P, W1, W2, W3, L = p

    M_plus_1, M_plus_2, M_plus_3, M_plus = SDIComposite.C2022.Appendix2_simple_span_Eq_C_A2_1_to_C_A2_3(P, W1, W2, W3, L)

    return M_plus_1, M_plus_2, M_plus_3, M_plus

end

function double_span_M(p)

    P, W1, W2, W3, L = p

    M_plus_1, M_plus_2, M_plus_3, M_plus, M_neg_1, M_neg_2, M_neg = SDIComposite.C2022.Appendix2_double_span_Eq_C_A2_8_to_C_A2_12(P, W1, W2, W3, L)

    return M_plus_1, M_plus_2, M_plus_3, M_plus, M_neg_1, M_neg_2, M_neg

end


function triple_span_M(p)

    P, W1, W2, W3, L = p

    M_plus_1, M_plus_2, M_plus_3, M_plus, M_neg_1, M_neg_2, M_neg = SDIComposite.C2022.Appendix2_triple_span_Eq_C_A2_20_to_C_A2_24(P, W1, W2, W3, L)

    return M_plus_1, M_plus_2, M_plus_3, M_plus, M_neg_1, M_neg_2, M_neg

end


function simple_span_web_crippling(p)

    P, W1, W2, W3, L = p

    P_ext_1, P_ext_2, P_ext_3 = SDIComposite.C2022.Appendix2_simple_span_Eq_C_A2_4_to_C_A2_6(W1, W2, W3, L, P)

    return P_ext_1, P_ext_2, P_ext_3

end

function double_span_web_crippling(p)

    P, W1, W2, W3, L = p

    P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3 = SDIComposite.C2022.Appendix2_double_span_Eq_C_A2_13_to_C_A2_18(W1, W2, W3, L, P)

    return P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3

end


function triple_span_web_crippling(p)

    P, W1, W2, W3, L = p

    P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3 = SDIComposite.C2022.Appendix2_triple_span_Eq_C_A2_25_to_C_A2_30(W1, W2, W3, L, P)

    return P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3

end


function simple_span_Δ(p)

    W1, L, E, I = p

    Δ = SDIComposite.C2022.Appendix2_simple_span_Eq_C_A2__7(W1, L, E, I)

    return Δ

end


function double_span_Δ(p)

    W1, L, E, I = p

    Δ = SDIComposite.C2022.Appendix2_double_span_Eq_C_A2_19(W1, L, E, I)

    return Δ

end


function triple_span_Δ(p)

    W1, L, E, I = p

    Δ = SDIComposite.C2022.Appendix2_triple_span_Eq_C_A2_31(W1, L, E, I)

    return Δ

end


function simple_span_M_V(p)

    P, W1, W2, W3, aMn, aVn, L = p

    M_plus_1, M_plus_2, M_plus_3, M_plus = SDIComposite.C2022.Appendix2_simple_span_Eq_C_A2_1_to_C_A2_3(P, W1, W2, W3, L)

    P_ext_1, P_ext_2, P_ext_3 = SDIComposite.C2022.Appendix2_simple_span_Eq_C_A2_4_to_C_A2_6(W1, W2, W3, L, P)

    Mbar = M_plus 
    Vbar = max(P_ext_1, P_ext_2, P_ext_3)  #this is W1+W2 loading condition, divide by 2 to convert reaction to shear 

    Maℓo = aMn 
    Va = aVn 
    
    interaction = AISIS100.v16S3.h21(Mbar, Vbar, Maℓo, Va)

    return interaction

end



function double_span_M_V(p)

    P, W1, W2, W3, aMn, aVn, L = p

    M_plus_1, M_plus_2, M_plus_3, M_plus, M_neg_1, M_neg_2, M_neg = SDIComposite.C2022.Appendix2_double_span_Eq_C_A2_8_to_C_A2_12(P, W1, W2, W3, L)

    P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3 = SDIComposite.C2022.Appendix2_double_span_Eq_C_A2_13_to_C_A2_18(W1, W2, W3, L, P)

    Mbar = max(M_neg, M_plus) 
    Vbar = max(P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3)  #this is the maximum reaction, which occurs at the support with the largest shear demand 

    Maℓo = aMn 
    Va = aVn 
    
    interaction = AISIS100.v16S3.h21(Mbar, Vbar, Maℓo, Va)

    return interaction

end


function triple_span_M_V(p)

    P, W1, W2, W3, aMn, aVn, L = p

    M_plus_1, M_plus_2, M_plus_3, M_plus, M_neg_1, M_neg_2, M_neg = SDIComposite.C2022.Appendix2_triple_span_Eq_C_A2_20_to_C_A2_24(P, W1, W2, W3, L)

    P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3 = SDIComposite.C2022.Appendix2_triple_span_Eq_C_A2_25_to_C_A2_30(W1, W2, W3, L, P)

    # https://faculty-legacy.arch.tamu.edu/anichols/index_files/courses/arch331/NS8-2beamdiagrams.pdf
    Mbar = max(M_neg, M_plus)
    Vbar = max(P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3)  #this is the maximum reaction, which occurs at the support with the largest shear demand

    Maℓo = aMn
    Va = aVn

    interaction = AISIS100.v16S3.h21(Mbar, Vbar, Maℓo, Va)

    return interaction

end


function construction_loads(inputs::ConstructionLoadInputs, configuration::AbstractString)

    p_load       = (inputs.P, inputs.W1, inputs.W2, inputs.W3, inputs.L)
    p_deflection = (inputs.W1, inputs.L, inputs.E, inputs.I)

    if configuration == "Single"

        M_plus_1, M_plus_2, M_plus_3, M_plus = simple_span_M(p_load)
        M_neg_1 = nothing; M_neg_2 = nothing; M_neg = nothing

        P_ext_1, P_ext_2, P_ext_3 = simple_span_web_crippling(p_load)
        P_int_1 = nothing; P_int_2 = nothing; P_int_3 = nothing

        Δ = simple_span_Δ(p_deflection)

        Mbar = M_plus
        Vbar = max(P_ext_1, P_ext_2, P_ext_3)

    elseif configuration == "Double"

        M_plus_1, M_plus_2, M_plus_3, M_plus, M_neg_1, M_neg_2, M_neg = double_span_M(p_load)

        P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3 = double_span_web_crippling(p_load)

        Δ = double_span_Δ(p_deflection)

        Mbar = max(M_neg, M_plus)
        Vbar = max(P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3)

    elseif configuration == "Triple"

        M_plus_1, M_plus_2, M_plus_3, M_plus, M_neg_1, M_neg_2, M_neg = triple_span_M(p_load)

        P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3 = triple_span_web_crippling(p_load)

        Δ = triple_span_Δ(p_deflection)

        Mbar = max(M_neg, M_plus)
        Vbar = max(P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3)

    else
        error("Configuration must be \"Single\", \"Double\", or \"Triple\"")
    end

    interaction = AISIS100.v16S3.h21(Mbar / inputs.unit_width, Vbar / inputs.unit_width, inputs.aMn, inputs.aVn)

    return ConstructionLoadOutputs(
        inputs,
        M_plus_1, M_plus_2, M_plus_3, M_plus,
        M_neg_1, M_neg_2, M_neg,
        P_ext_1, P_ext_2, P_ext_3,
        P_int_1, P_int_2, P_int_3,
        Δ,
        interaction
    )

end
