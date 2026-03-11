
struct BareShearInputs

    t
    Fy
    E
    μ
    kv

    steel_deck_depth
    steel_deck_web_angle_from_horizon
    steel_deck_outside_radius
    number_of_webs_per_panel

    unit_width

    design_method

end



struct BareShearOutputs

    inputs

    h

    Aw
    Vy

    Fcr
    Vcr

    Vn_web

    aVn_web

    Vn_unit

    aVn_unit

end



function calculate_bare_shear_strength(inputs, S100_version)

    (;     
    t, 
    Fy, 
    E, 
    μ,
    kv, 

    steel_deck_depth, 
    steel_deck_web_angle_from_horizon,
    steel_deck_outside_radius, 
    number_of_webs_per_panel,

    unit_width,

    design_method,

    ) = inputs


    h = steel_deck_depth / sind(steel_deck_web_angle_from_horizon) - 2 * steel_deck_outside_radius


    Aw, Vy = AISIS100.v16S3.g215_6(h, t, Fy)

    Fcr = AISIS100.v16S3.g232(E, μ, kv, h, t)

    Vcr = AISIS100.v16S3.g231(h, t, Fcr)

    if S100_version == "v24"
        Vn_web, aVn_web = v2024.g2_1__1_2_3(Vcr, Vy, design_method)
    else  # v16
        Vn_web, aVn_web = AISIS100.v16S3.g2_1__1_2_3(Vcr, Vy, "AISI S100-16 $design_method")
    end

    Vn_unit  = Vn_web  * number_of_webs_per_panel / unit_width
    aVn_unit = aVn_web * number_of_webs_per_panel / unit_width

    outputs = BareShearOutputs(
        inputs,

        h,

        Aw,
        Vy,

        Fcr,
        Vcr,

        Vn_web,

        aVn_web,

        Vn_unit,

        aVn_unit,

    )


    return outputs 

end