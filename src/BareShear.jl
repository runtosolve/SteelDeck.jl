
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

end



struct BareShearOutputs
        
    inputs 

    h 

    Aw 
    Vy 

    Fcr
    Vcr 

    Vn_web 
    aVn_web_ASD 
    aVn_web_LRFD 

    Vn_unit 
    aVn_unit_ASD
    aVn_unit_LRFD 

end



function calculate_bare_shear_strength(inputs)

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

    unit_width 

    ) = inputs 


    h = steel_deck_depth / sind(steel_deck_web_angle_from_horizon) - 2 * steel_deck_outside_radius


    Aw, Vy = AISIS100.v16S3.g215_6(h, t, Fy)

    Fcr = AISIS100.v16S3.g232(E, μ, kv, h, t)

    Vcr = AISIS100.v16S3.g231(h, t, Fcr)

    design_code = "AISI S100-16 ASD"
    Vn_web, aVn_web_ASD = AISIS100.v16S3.g2_1__1_2_3(Vcr, Vy, design_code)

    design_code = "AISI S100-16 LRFD"
    Vn_web, aVn_web_LRFD = AISIS100.v16S3.g2_1__1_2_3(Vcr, Vy, design_code)

    Vn_unit = Vn_web * number_of_webs_per_panel / unit_width
    aVn_unit_ASD = aVn_web_ASD * number_of_webs_per_panel / unit_width
    aVn_unit_LRFD = aVn_web_LRFD * number_of_webs_per_panel / unit_width


    outputs =    BareShearOutputs(
        inputs, 

        h, 

        Aw, 
        Vy, 

        Fcr,
        Vcr, 

        Vn_web, 
        aVn_web_ASD, 
        aVn_web_LRFD, 

        Vn_unit, 
        aVn_unit_ASD,
        aVn_unit_LRFD 
    )


    return outputs 

end