struct OneWayShearInputs

    t
    steel_deck_depth
    steel_deck_web_angle_from_horizon
    steel_deck_radius
    steel_deck_trough_width
    number_of_troughs_in_panel
    unit_width

    total_slab_depth

    Fy
    E
    μ

    kv

    λ
    fc

    design_method

end



struct OneWayShearOutputs

    inputs

    h 
    Aw
    Vy
    Fcr
    Vcr
    VD
    
    Δ_pitch 
    Ac 
    Vc 
    aVn_per_pitch
    aVn_unit 

end




function calculate_one_way_composite_shear_strength(inputs)

    (; t, 
    steel_deck_depth,
    steel_deck_web_angle_from_horizon,
    steel_deck_radius,
    steel_deck_trough_width, 
    number_of_troughs_in_panel,
    unit_width,
    
    total_slab_depth, 

    Fy,
    E, 
    μ,

    kv,

    λ,
    fc,
    design_method) = inputs

    ####steel strength 

    h = steel_deck_depth / sind(steel_deck_web_angle_from_horizon) - 2 * steel_deck_radius


    Aw, Vy = AISIS100.v16S3.g215_6(h, t, Fy)

    Fcr = AISIS100.v16S3.g232(E, μ, kv, h, t)

    Vcr = AISIS100.v16S3.g231(h, t, Fcr)

    design_code = "nominal"
    VD, aVD = AISIS100.v16S3.g2_1__1_2_3(Vcr, Vy, design_code)


    #####concrete strength 

    Δ_pitch = total_slab_depth * tand(90.0 - steel_deck_web_angle_from_horizon)
    Ac = mean([steel_deck_trough_width, steel_deck_trough_width + 2 * Δ_pitch]) * total_slab_depth


    Vc = SDIComposite.C2017.Eq2_4_8_a(λ, fc, Ac)

    if design_method == "ASD"
        aVn_per_pitch = SDIComposite.C2017.Eq2_4_7c(Vc, 2 * VD, fc, Ac) #2 * VD for two webs per pitch 
    elseif design_method == "LRFD"
        aVn_per_pitch = SDIComposite.C2017.Eq2_4_7a(Vc, 2 * VD, fc, Ac) #2 * VD for two webs per pitch
    end

    aVn_unit = (aVn_per_pitch * number_of_troughs_in_panel) / unit_width

    #package outputs 

    outputs = SteelDeck.OneWayShearOutputs(

                inputs,

                h, 
                Aw,
                Vy,
                Fcr,
                Vcr,
                VD,

                Δ_pitch, 
                Ac,
                Vc,
                aVn_per_pitch,
                aVn_unit 
    )

    return outputs 

end