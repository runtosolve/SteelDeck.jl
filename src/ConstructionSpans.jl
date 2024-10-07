struct ConstructionSpanInputs

    bare_properties 
    panel_width
    unit_width 
    Dw
    ph
    deck_flute_widths 
    dd 
    number_empty_flutes 
    γs
    Es
    Fy 

    γc
    h 

    wlc 
    Plc 

    Δ_over_L_limit

    yc_steel 
    Ixx_steel 
    As 
    Wc 
    fc 
    b 
    β1 

end



struct ConstructionSpanOutputs 

    Ec
    n 
    hc 
    d 
    ρ
    ycc 
    ycs 
    Isf 
    Icr 
    K1 
    K3 
    
    flexural_properties 

    flute_area 
    W1_concrete
    W1_steel 
    W1 

    W2 
    P 

    aMn 

    ℓ_M_plus

    E 
    I 

    ℓ_Δ

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

function triple_span_M_plus(u, p)

    P, W1, W2, aMn = p

    M_plus_1, M_plus_2, M_plus, M_neg = SDIComposite.C2017.Appendix1_Figure_1_triple_span(P, W1, W2, u[1])

    return (M_plus - aMn)

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

function calculate_M_plus_spans(P, W1, W2, aMn)

    ℓ_M_plus = OrderedDict{String,Float64}()

    uspan = [0.00001, 99999.0]
    p = (P, W1, W2, aMn)

    prob = IntervalNonlinearProblem(simple_span_M_plus, uspan, p)
    sol = solve(prob)
    ℓ_M_plus["simple span"] = sol.u

    prob = IntervalNonlinearProblem(double_span_M_plus, uspan, p)
    sol = solve(prob)
    ℓ_M_plus["double span"] = sol.u

    prob = IntervalNonlinearProblem(triple_span_M_plus, uspan, p)
    sol = solve(prob)
    ℓ_M_plus["triple span"] = sol.u

    return ℓ_M_plus

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


function calculate_construction_spans(inputs)


    (;   bare_properties, 
    panel_width,
    unit_width, 
    Dw,
    ph,
    deck_flute_widths, 
    dd, 
    number_empty_flutes, 
    γs,
    Es, 

    γc,
    h, 

    wlc, 
    Plc, 

    Δ_over_L_limit,

    yc_steel, 
    Ixx_steel, 
    As,
    Wc, 
    fc, 
    b, 
    β1
    ) = inputs 


    Ec = Wc^1.5 * (fc/1000)^0.5 * 1000

    n = Es / Ec

    hc = h - dd 

    unit_ratio = unit_width / panel_width

    d = h - yc_steel

    ρ = As / (b * d)

    ycc = SDIComposite.C2017.EqA5__1(d, ρ, n, hc)

    ycs = d - ycc 

    Isf = Ixx_steel * unit_ratio

    Icr = SDIComposite.C2017.Eq_A5__2(b, n, ycc, As, ycs, Isf)


    K1 = SDIComposite.C2017.EqA2__12(Dw, ph)
    K3 = 1.4


    ####flexural strength 
    inputs = (Es, Fy, fc, β1, d, dd, b, h, ρ, As, Icr, ycc, K1, K3)
    flexural_properties = SDIComposite.C2017.SectionA2__2_flexural_strength(inputs)

    flute_area = mean(deck_flute_widths) * dd
    W1_concrete = (h * panel_width - flute_area * number_empty_flutes) / 144 * γc / (panel_width / 12)
    W1_steel = As / unit_width / 12 * γs
    W1 = W1_concrete + W1_steel

    W2 = wlc 
    P = Plc 


    if isnothing(flexural_properties.Mno_ASD)
        aMn = flexural_properties.Mro_ASD
    else
        aMn = flexural_properties.Mno_ASD
    end

    #convert from lbf-in to lbf-ft
    aMn = aMn / 12

    ###########


    ℓ_M_plus = SteelDeck.calculate_M_plus_spans(P, W1, W2, aMn)


    E = Es * 12.0^2  #convert to psf 
    I = bare_properties.Ixx_unit / 12^4 #convert to ft^4

    ℓ_Δ = SteelDeck.calculate_Δ_spans(W1, E, I, Δ_over_L_limit)

    outputs = ConstructionSpanOutputs(
                                    Ec,
                                    n,
                                    hc, 
                                    d, 
                                    ρ,
                                    ycc, 
                                    ycs, 
                                    Isf, 
                                    Icr, 
                                    K1, 
                                    K3, 
                                    
                                    flexural_properties, 

                                    flute_area, 
                                    W1_concrete,
                                    W1_steel, 
                                    W1, 

                                    W2, 
                                    P, 

                                    aMn, 

                                    ℓ_M_plus,

                                    E, 
                                    I, 

                                    ℓ_Δ)


    return outputs 

end