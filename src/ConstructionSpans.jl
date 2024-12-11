
struct ConstructionSpanInputs 

    bare_deck_properties 
    web_crippling_properties
    bare_shear_properties 

    E 
    deck_flute_widths 
    number_empty_flutes 
    dd 
    panel_width 
    γs

    γc

    h 

    unit_width 

    wlc 
    Plc 
    wcdl

    Δ_over_L_limit

end 


struct ConstructionSpanOutputs 

    inputs

    flute_area
    As 
    I  

    W1_concrete
    W1_steel
    
    W1_lc 
    W2_lc 
    P_lc 
    
    W1_cdl 
    W2_cdl
    P_cdl  

    aMnℓ_pos_unit_ASD
    ℓ_M_plus_lc
    ℓ_M_plus_cdl

    aMnℓ_neg_unit_ASD
    ℓ_M_neg_lc
    ℓ_M_neg_cdl

    aPn_OFE
    aPn_TFE
    aPn_E 
    ℓ_web_crippling_ext_lc
    ℓ_web_crippling_ext_cdl

    aPn_OFI
    aPn_TFI
    aPn_I
    ℓ_web_crippling_int_lc
    ℓ_web_crippling_int_cdl

    ℓ_Δ_lc
    ℓ_Δ_cdl

    aVn_unit_ASD
    ℓ_M_V_lc
    ℓ_M_V_cdl

    ℓ_min_lc
    ℓ_min_cdl
    controlling_limit_states_lc
    controlling_limit_states_cdl

    ℓ_min
    controlling_load_combination
 
end



function simple_span_M_plus(u, p)

    P, W1, W2, aMn = p

    M_plus_1, M_plus_2, M_plus = SDIComposite.C2017.Appendix1_Figure_1_simple_span(P, W1, W2, u[1])

    return (M_plus - aMn)

end

function double_span_M_plus(u, p)

    P, W1, W2, aMn = p

    M_plus_1, M_plus_2, M_plus, M_neg = SDIComposite.C2017.Appendix1_Figure_1_double_span(P, W1, W2, u[1])

    return (M_plus - aMn)

end

function double_span_M_neg(u, p)

    P, W1, W2, aMn = p

    M_plus_1, M_plus_2, M_plus, M_neg = SDIComposite.C2017.Appendix1_Figure_1_double_span(P, W1, W2, u[1])

    return (M_neg - aMn)

end

function triple_span_M_plus(u, p)

    P, W1, W2, aMn = p

    M_plus_1, M_plus_2, M_plus, M_neg = SDIComposite.C2017.Appendix1_Figure_1_triple_span(P, W1, W2, u[1])

    return (M_plus - aMn)

end

function triple_span_M_neg(u, p)

    P, W1, W2, aMn = p

    M_plus_1, M_plus_2, M_plus, M_neg = SDIComposite.C2017.Appendix1_Figure_1_triple_span(P, W1, W2, u[1])

    return (M_neg - aMn)

end


function simple_span_web_crippling_ext(u, p)

    P, W1, W2, aPn = p

    P_ext_1, P_ext_2 = SDIComposite.C2017.Appendix1_Figure_2_simple_span(W1, W2, u[1], P)

    return (maximum([P_ext_1, P_ext_2]) - aPn)

end

function double_span_web_crippling_ext(u, p)

    P, W1, W2, aPn = p

    P_ext_1, P_ext_2, P_int_1, P_int_2 = SDIComposite.C2017.Appendix1_Figure_2_double_span(W1, W2, u[1], P)

    return (maximum([P_ext_1, P_ext_2]) - aPn)

end

function double_span_web_crippling_int(u, p)

    P, W1, W2, aPn = p

    P_ext_1, P_ext_2, P_int_1, P_int_2 = SDIComposite.C2017.Appendix1_Figure_2_double_span(W1, W2, u[1], P)

    return (maximum([P_int_1, P_int_2]) - aPn)

end



function triple_span_web_crippling_ext(u, p)

    P, W1, W2, aPn = p

    P_ext_1, P_ext_2, P_int_1, P_int_2 = SDIComposite.C2017.Appendix1_Figure_2_triple_span(W1, W2, u[1], P)

    return (maximum([P_ext_1, P_ext_2]) - aPn)

end


function triple_span_web_crippling_int(u, p)

    P, W1, W2, aPn = p

    P_ext_1, P_ext_2, P_int_1, P_int_2 = SDIComposite.C2017.Appendix1_Figure_2_triple_span(W1, W2, u[1], P)

    return (maximum([P_int_1, P_int_2]) - aPn)

end


function simple_span_Δ(u, p)

    W1, E, I, Δ_over_L_limit = p

    Δ = SDIComposite.C2017.Appendix1_Figure_3_simple_span(W1, u[1], E, I)

    return (Δ - Δ_over_L_limit * u[1])

end


function double_span_Δ(u, p)

    W1, E, I, Δ_over_L_limit = p

    Δ = SDIComposite.C2017.Appendix1_Figure_3_double_span(W1, u[1], E, I)

    return (Δ - Δ_over_L_limit * u[1])

end


function triple_span_Δ(u, p)

    W1, E, I, Δ_over_L_limit = p

    Δ = SDIComposite.C2017.Appendix1_Figure_3_triple_span(W1, u[1], E, I)

    return (Δ - Δ_over_L_limit * u[1])

end


function double_span_M_V(u, p)

    P, W1, W2, aMn, aVn = p

    M_plus_1, M_plus_2, M_plus, M_neg = SDIComposite.C2017.Appendix1_Figure_1_double_span(P, W1, W2, u[1])

    # P_ext_1, P_ext_2, P_int_1, P_int_2 = SDIComposite.C2017.Appendix1_Figure_2_double_span(W1, W2, u[1], P)

    Mbar = M_neg 
    Vbar = 0.50 * (W1 + W2) * u[1]  #this is W1+W2 loading condition, divide by 2 to convert reaction to shear 

    Maℓo = aMn 
    Va = aVn 
    
    interaction = AISIS100.v16S3.h21(Mbar, Vbar, Maℓo, Va)

    return (interaction - 1.0)

end


function triple_span_M_V(u, p)

    P, W1, W2, aMn, aVn = p

    M_plus_1, M_plus_2, M_plus, M_neg = SDIComposite.C2017.Appendix1_Figure_1_triple_span(P, W1, W2, u[1])

    # P_ext_1, P_ext_2, P_int_1, P_int_2 = SDIComposite.C2017.Appendix1_Figure_2_triple_span(W1, W2, u[1], P)

    # https://faculty-legacy.arch.tamu.edu/anichols/index_files/courses/arch331/NS8-2beamdiagrams.pdf
    Mbar = M_neg 
    Vbar = 0.617 * (W1 + W2) * u[1]  #this is W1+W2 loading condition 

    Maℓo = aMn 
    Va = aVn 
    
    interaction = AISIS100.v16S3.h21(Mbar, Vbar, Maℓo, Va)

    return (interaction - 1.0)

end





function calculate_M_plus_spans(P, W1, W2, aMn)

    ℓ = OrderedDict{String,Float64}()

    uspan = [0.00001, 99999.0]
    p = (P, W1, W2, aMn)

    prob = IntervalNonlinearProblem(simple_span_M_plus, uspan, p)
    sol = solve(prob)
    ℓ["simple span"] = sol.u

    prob = IntervalNonlinearProblem(double_span_M_plus, uspan, p)
    sol = solve(prob)
    ℓ["double span"] = sol.u

    prob = IntervalNonlinearProblem(triple_span_M_plus, uspan, p)
    sol = solve(prob)
    ℓ["triple span"] = sol.u

    return ℓ

end

function calculate_M_neg_spans(P, W1, W2, aMn)

    ℓ = OrderedDict{String,Float64}()

    uspan = [0.00001, 99999.0]
    p = (P, W1, W2, aMn)

    prob = IntervalNonlinearProblem(double_span_M_neg, uspan, p)
    sol = solve(prob)
    ℓ["double span"] = sol.u

    prob = IntervalNonlinearProblem(triple_span_M_neg, uspan, p)
    sol = solve(prob)
    ℓ["triple span"] = sol.u

    return ℓ

end


function calculate_web_crippling_ext_spans(P, W1, W2, aPn)

    ℓ = OrderedDict{String,Float64}()

    uspan = [0.00001, 99999.0]
    p = (P, W1, W2, aPn)

    prob = IntervalNonlinearProblem(simple_span_web_crippling_ext, uspan, p)
    sol = solve(prob)
    ℓ["simple span"] = sol.u

    prob = IntervalNonlinearProblem(double_span_web_crippling_ext, uspan, p)
    sol = solve(prob)
    ℓ["double span"] = sol.u

    prob = IntervalNonlinearProblem(triple_span_web_crippling_ext, uspan, p)
    sol = solve(prob)
    ℓ["triple span"] = sol.u

    return ℓ

end


function calculate_web_crippling_int_spans(P, W1, W2, aPn)

    ℓ = OrderedDict{String,Float64}()

    uspan = [0.00001, 99999.0]
    p = (P, W1, W2, aPn)


    prob = IntervalNonlinearProblem(double_span_web_crippling_int, uspan, p)
    sol = solve(prob)
    ℓ["double span"] = sol.u

    prob = IntervalNonlinearProblem(triple_span_web_crippling_int, uspan, p)
    sol = solve(prob)
    ℓ["triple span"] = sol.u

    return ℓ

end



function calculate_Δ_spans(W1, E, I, Δ_over_L_limit)

    ℓ_Δ = OrderedDict{String,Float64}()

    uspan = [0.00001, 99999.0]
    p = (W1, E, I, Δ_over_L_limit)
    prob = IntervalNonlinearProblem(simple_span_Δ, uspan, p)
    sol = solve(prob)
    ℓ_Δ["simple span"] = sol.u

    prob = IntervalNonlinearProblem(double_span_Δ, uspan, p)
    sol = solve(prob)
    ℓ_Δ["double span"] = sol.u

    prob = IntervalNonlinearProblem(triple_span_Δ, uspan, p)
    sol = solve(prob)
    ℓ_Δ["triple span"] = sol.u

    return ℓ_Δ

end


function calculate_M_V_spans(P, W1, W2, aMn, aVn)

    ℓ = OrderedDict{String,Float64}()

    uspan = [0.00001, 99999.0]
    p = (P, W1, W2, aMn, aVn)

    prob = IntervalNonlinearProblem(double_span_M_V, uspan, p)
    sol = solve(prob)
    ℓ["double span"] = sol.u


    prob = IntervalNonlinearProblem(triple_span_M_V, uspan, p)
    sol = solve(prob)
    ℓ["triple span"] = sol.u

    return ℓ

end




function calculate_construction_spans(inputs)

    (;
        bare_deck_properties, 
        web_crippling_properties,
        bare_shear_properties, 
    
        E, 
        deck_flute_widths, 
        number_empty_flutes, 
        dd, 
        panel_width, 
        γs,
    
        γc,
    
        h, 
    
        unit_width, 
    
        wlc, 
        Plc, 
        wcdl,
    
        Δ_over_L_limit


    ) = inputs



    As_whole_panel = bare_deck_properties.section_properties.A
    Ixx_steel = bare_deck_properties.section_properties.Ixx

    flute_area = mean(deck_flute_widths) * dd
    W1_concrete = (h * panel_width - flute_area * number_empty_flutes) / 144 * γc / (panel_width / 12)
    As = As_whole_panel * (unit_width / panel_width)
    I = Ixx_steel * (unit_width / panel_width)
    W1_steel = As / unit_width / 12 * γs
    
    ####  SDI C-2017, Eq. 2.4.1, Eq. 2.4.2 considered 
    
    W1 = W1_concrete + W1_steel
    W2 = wlc 
    P = Plc 

    aMnℓ_pos_unit_ASD = bare_deck_properties.aMnℓ_pos_unit_ASD
    ℓ_M_plus_lc = calculate_M_plus_spans(P, W1/12, W2/12, aMnℓ_pos_unit_ASD)

    aMnℓ_neg_unit_ASD = bare_deck_properties.aMnℓ_neg_unit_ASD
    ℓ_M_neg_lc = calculate_M_neg_spans(P, W1/12, W2/12, aMnℓ_neg_unit_ASD)

    aPn_OFE = web_crippling_properties.OFE.aPn_unit
    aPn_TFE = web_crippling_properties.TFE.aPn_unit
    aPn_E = minimum([web_crippling_properties.OFE.aPn_unit, web_crippling_properties.TFE.aPn_unit])
    ℓ_web_crippling_ext_lc =  calculate_web_crippling_ext_spans(P, W1/12, W2/12, aPn_E)

    aPn_OFI = web_crippling_properties.OFI.aPn_unit
    aPn_TFI = web_crippling_properties.TFI.aPn_unit
    aPn_I = minimum([web_crippling_properties.OFI.aPn_unit, web_crippling_properties.TFI.aPn_unit])
    ℓ_web_crippling_int_lc =  calculate_web_crippling_int_spans(P, W1/12, W2/12, aPn_I)

    ℓ_Δ_lc = calculate_Δ_spans(W1/12, E, I, Δ_over_L_limit)

    aVn_unit_ASD = bare_shear_properties.aVn_unit_ASD
    aMnℓ_neg_unit_ASD = bare_deck_properties.aMnℓ_neg_unit_ASD
    ℓ_M_V_lc = calculate_M_V_spans(P, W1/12, W2/12, aMnℓ_neg_unit_ASD, aVn_unit_ASD)


    ℓ_min_lc, controlling_limit_states_lc = find_minimum_construction_spans(ℓ_M_plus_lc, ℓ_M_neg_lc, ℓ_web_crippling_ext_lc, ℓ_web_crippling_int_lc, ℓ_Δ_lc, ℓ_M_V_lc)


    W1_lc = W1
    W2_lc = W2 
    P_lc = P

    ##### SDI C-2017, Eq. 2.4.3 considered  

    W1 = W1_steel
    W2 = wcdl 
    P = 0.0 

    ℓ_M_plus_cdl = calculate_M_plus_spans(P, W1/12, W2/12, aMnℓ_pos_unit_ASD)

    ℓ_M_neg_cdl = calculate_M_neg_spans(P, W1/12, W2/12, aMnℓ_neg_unit_ASD)

    ℓ_web_crippling_ext_cdl =  calculate_web_crippling_ext_spans(P, W1/12, W2/12, aPn_E)

    ℓ_web_crippling_int_cdl =  calculate_web_crippling_int_spans(P, W1/12, W2/12, aPn_I)

    ℓ_Δ_cdl = calculate_Δ_spans(W1/12, E, I, Δ_over_L_limit)

    ℓ_M_V_cdl = calculate_M_V_spans(P, W1/12, W2/12, aMnℓ_neg_unit_ASD, aVn_unit_ASD)

    ℓ_min_cdl, controlling_limit_states_cdl = find_minimum_construction_spans(ℓ_M_plus_cdl, ℓ_M_neg_cdl, ℓ_web_crippling_ext_cdl, ℓ_web_crippling_int_cdl, ℓ_Δ_cdl, ℓ_M_V_cdl)

    W1_cdl = W1
    W2_cdl = W2 
    P_cdl = P


    ######

    load_combinations = ["Eq. 2.4.1, Eq. 2.4.2", "Eq. 2.4.3"]

    ℓ_min = OrderedDict{String, Float64}()
    controlling_load_combination = OrderedDict{String, String}()

    ℓ_min["simple span"] = minimum([ℓ_min_lc["simple span"], ℓ_min_cdl["simple span"]])
    ℓ_min["double span"] = minimum([ℓ_min_lc["double span"], ℓ_min_cdl["double span"]])
    ℓ_min["triple span"] = minimum([ℓ_min_lc["triple span"], ℓ_min_cdl["triple span"]])
    
    index = argmin([ℓ_min_lc["simple span"], ℓ_min_cdl["simple span"]])
    controlling_load_combination["simple span"] = load_combinations[index]

    index = argmin([ℓ_min_lc["double span"], ℓ_min_cdl["double span"]])
    controlling_load_combination["double span"] = load_combinations[index]

    index = argmin([ℓ_min_lc["triple span"], ℓ_min_cdl["triple span"]])
    controlling_load_combination["triple span"] = load_combinations[index]

 
    outputs = ConstructionSpanOutputs(

        inputs,

        flute_area,
        As, 
        I,  

        W1_concrete,
        W1_steel,
        
        W1_lc, 
        W2_lc, 
        P_lc, 

        W1_cdl, 
        W2_cdl, 
        P_cdl, 

        aMnℓ_pos_unit_ASD,
        ℓ_M_plus_lc,
        ℓ_M_plus_cdl,

        aMnℓ_neg_unit_ASD,
        ℓ_M_neg_lc,
        ℓ_M_neg_cdl,

        aPn_OFE,
        aPn_TFE,
        aPn_E, 
        ℓ_web_crippling_ext_lc,
        ℓ_web_crippling_ext_cdl,

        aPn_OFI,
        aPn_TFI,
        aPn_I,
        ℓ_web_crippling_int_lc,
        ℓ_web_crippling_int_cdl,

        ℓ_Δ_lc,
        ℓ_Δ_cdl,

        aVn_unit_ASD,
        ℓ_M_V_lc,
        ℓ_M_V_cdl,

        ℓ_min_lc,
        ℓ_min_cdl,
        controlling_limit_states_lc,
        controlling_limit_states_cdl,

        ℓ_min,
        controlling_load_combination
        


    )

    return outputs 

end


function find_minimum_construction_spans(ℓ_M_plus, ℓ_M_neg, ℓ_web_crippling_ext, ℓ_web_crippling_int, ℓ_Δ, ℓ_M_V)

    ℓ = OrderedDict{String, Float64}()

    controlling_limit_states = OrderedDict{String, String}()

    simple_span_limit_states = ["positive flexure", "web crippling", "deflection"]

    continuous_span_limit_states = ["positive flexure", "negative flexure", "web crippling exterior", "web crippling interior", "deflection", "negative flexure + shear"]

    ℓ["simple span"] = minimum([ℓ_M_plus["simple span"], ℓ_web_crippling_ext["simple span"], ℓ_Δ["simple span"]])
    index = argmin([ℓ_M_plus["simple span"], ℓ_web_crippling_ext["simple span"], ℓ_Δ["simple span"]])
    controlling_limit_states["simple span"] = simple_span_limit_states[index]

    ℓ["double span"] = minimum([ℓ_M_plus["double span"], ℓ_M_neg["double span"], ℓ_web_crippling_ext["double span"], ℓ_web_crippling_int["double span"], ℓ_Δ["double span"], ℓ_M_V["double span"]])
    index = argmin([ℓ_M_plus["double span"], ℓ_M_neg["double span"], ℓ_web_crippling_ext["double span"], ℓ_web_crippling_int["double span"], ℓ_Δ["double span"], ℓ_M_V["double span"]])
    controlling_limit_states["double span"] = continuous_span_limit_states[index]

    ℓ["triple span"] = minimum([ℓ_M_plus["triple span"], ℓ_M_neg["triple span"], ℓ_web_crippling_ext["triple span"], ℓ_web_crippling_int["triple span"], ℓ_Δ["triple span"], ℓ_M_V["triple span"]])
    index = argmin([ℓ_M_plus["triple span"], ℓ_M_neg["triple span"], ℓ_web_crippling_ext["triple span"], ℓ_web_crippling_int["triple span"], ℓ_Δ["triple span"], ℓ_M_V["triple span"]])
    controlling_limit_states["triple span"] = continuous_span_limit_states[index]

    return ℓ, controlling_limit_states

end