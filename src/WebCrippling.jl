struct WebCripplingInputs

    t
    steel_deck_depth
    steel_deck_web_angle_from_horizon
    steel_deck_outside_radius 
    number_of_webs_per_panel

    unit_width 

    bearing_width 
    N
  
    Fy

    C 
    C_R 
    C_N
    C_h

    Ω_w
    ϕ_w
    ϕ_w_LSD
    
    design_code

end


struct WebCripplingOutputs

    inputs 

    h 
    θ
    R 

    Pn 
    aPn 
    aPn_unit

end

struct AllWebCripplingConditions

            OFE
            TFE
            OFI
            TFI

end




function calculate_web_crippling_strength(inputs)

    (; t,
    steel_deck_depth,
    steel_deck_web_angle_from_horizon,
    steel_deck_outside_radius, 
    number_of_webs_per_panel,

    unit_width, 
    
    bearing_width,
    N,
  
    Fy,

    C, 
    C_R, 
    C_N,
    C_h,

    Ω_w,
    ϕ_w,
    ϕ_w_LSD,
    
    design_code
    ) = inputs 


    h = steel_deck_depth / sind(steel_deck_web_angle_from_horizon) - 2 * steel_deck_outside_radius
    θ = steel_deck_web_angle_from_horizon
    R = steel_deck_outside_radius - t

    Pn, aPn = AISIS100.v16S3.g51(t, h, Fy, θ, C, C_R, R, C_N, N, C_h, ϕ_w, Ω_w, ϕ_w_LSD, design_code)

    aPn_unit = aPn * number_of_webs_per_panel / unit_width

    


    outputs = WebCripplingOutputs(

                inputs,

                h, 
                θ,
                R, 

                Pn, 
                aPn,
                aPn_unit 
    )

    return outputs 

end






function calculate_all_web_crippling_conditions(inputs)

    (; t,
    steel_deck_depth,
    steel_deck_web_angle_from_horizon,
    steel_deck_outside_radius, 
    number_of_webs_per_panel,

    unit_width, 
    
    
    
    bearing_width,
    N,
  
    Fy,

    C, 
    C_R, 
    C_N,
    C_h,

    Ω_w,
    ϕ_w,
    ϕ_w_LSD,
    
    design_code
    ) = inputs  



    #AISI S100-16 Table G5-5

    #One-Flange End

    N = bearing_width / 2

    C = 4.0
    C_R = 0.04
    C_N = 0.25
    C_h = 0.025

    Ω_w = 1.70
    ϕ_w = 0.90
    ϕ_w_LSD = 0.80 


    inputs = WebCripplingInputs(

        t,
        steel_deck_depth,
        steel_deck_web_angle_from_horizon,
        steel_deck_outside_radius,
        number_of_webs_per_panel,

        unit_width, 
         
        
        bearing_width,
        N, 
    
        Fy,

        C, 
        C_R,
        C_N,
        C_h,

        Ω_w,
        ϕ_w,
        ϕ_w_LSD,
        
        design_code,

    )


    OFE = calculate_web_crippling_strength(inputs)




    #Two-Flange End

    N = bearing_width / 2 

    C = 9.0
    C_R = 0.12
    C_N = 0.14
    C_h = 0.040

    Ω_w = 1.80
    ϕ_w = 0.85
    ϕ_w_LSD = 0.70 



    inputs = WebCripplingInputs(

        t,
        steel_deck_depth,
        steel_deck_web_angle_from_horizon,
        steel_deck_outside_radius,
        number_of_webs_per_panel,

        unit_width, 
         
        
        bearing_width,
        N, 
    
        Fy,

        C, 
        C_R,
        C_N,
        C_h,

        Ω_w,
        ϕ_w,
        ϕ_w_LSD,
        
        design_code,

    )

    TFE = calculate_web_crippling_strength(inputs)



    #One-Flange Interior

    N = bearing_width 

    C = 8.0
    C_R = 0.10
    C_N = 0.17
    C_h = 0.004

    Ω_w = 1.75
    ϕ_w = 0.85
    ϕ_w_LSD = 0.75 



    inputs = WebCripplingInputs(

        t,
        steel_deck_depth,
        steel_deck_web_angle_from_horizon,
        steel_deck_outside_radius, 
        number_of_webs_per_panel,

        unit_width, 
        
        
        bearing_width, 
        N,
    
        Fy,

        C, 
        C_R,
        C_N,
        C_h,

        Ω_w,
        ϕ_w,
        ϕ_w_LSD,
        
        design_code,

    )

    OFI = calculate_web_crippling_strength(inputs)


    #Two-Flange Interior

    N = bearing_width 

    C = 10.0
    C_R = 0.11
    C_N = 0.21
    C_h = 0.02

    Ω_w = 1.75
    ϕ_w = 0.85
    ϕ_w_LSD = 0.75 



    inputs = WebCripplingInputs(

        t,
        steel_deck_depth,
        steel_deck_web_angle_from_horizon,
        steel_deck_outside_radius,
        number_of_webs_per_panel,

        unit_width, 
         
        
        bearing_width, 
        N,
    
        Fy,

        C, 
        C_R,
        C_N,
        C_h,

        Ω_w,
        ϕ_w,
        ϕ_w_LSD,
        
        design_code,

    )

    TFI = calculate_web_crippling_strength(inputs)


    outputs = AllWebCripplingConditions(

            OFE,
            TFE,
            OFI,
            TFI

    )

    return outputs 

end
