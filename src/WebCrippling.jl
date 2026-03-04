struct WebCripplingInputs

    t
    fy
    panel_depth
    steel_deck_web_angle_from_horizon
    steel_deck_outside_radius
    number_of_webs_per_panel

    unit_width

    bearing_width

end


struct WebCripplingConditionInputs

    t
    fy
    panel_depth
    steel_deck_web_angle_from_horizon
    steel_deck_outside_radius
    number_of_webs_per_panel

    unit_width

    bearing_width

    N
    C
    C_R
    C_N
    C_h
    Ω_w
    ϕ_w
    ϕ_w_LSD

end


struct WebCripplingOutputs

    inputs

    h
    θ
    R

    Pn
    aPn_ASD
    aPn_LRFD
    aPn_unit_ASD
    aPn_unit_LRFD
end

struct AllWebCripplingConditions

            OFE
            TFE
            OFI
            TFI

end




function calculate_web_crippling_strength(inputs)

    (; t,
    fy,
    panel_depth,
    steel_deck_web_angle_from_horizon,
    steel_deck_outside_radius,
    number_of_webs_per_panel,

    unit_width,

    bearing_width,

    N,
    C,
    C_R,
    C_N,
    C_h,
    Ω_w,
    ϕ_w,
    ϕ_w_LSD,

    ) = inputs


    h = panel_depth / sind(steel_deck_web_angle_from_horizon) - 2 * steel_deck_outside_radius
    θ = steel_deck_web_angle_from_horizon
    R = steel_deck_outside_radius - t


    design_code = "ASD"
    Pn, aPn_ASD = AISIS100.v2024.g51(t, h, fy, θ, C, C_R, R, C_N, N, C_h, ϕ_w, Ω_w, ϕ_w_LSD, design_code)
    aPn_unit_ASD = aPn_ASD * number_of_webs_per_panel / unit_width

    design_code = "LRFD"
    Pn, aPn_LRFD = AISIS100.v2024.g51(t, h, fy, θ, C, C_R, R, C_N, N, C_h, ϕ_w, Ω_w, ϕ_w_LSD, design_code)
    aPn_unit_LRFD = aPn_LRFD * number_of_webs_per_panel / unit_width


    outputs = WebCripplingOutputs(

                inputs,

                h,
                θ,
                R,

                Pn,
                aPn_ASD,
                aPn_LRFD,
                aPn_unit_ASD,
                aPn_unit_LRFD,
    )

    return outputs

end




function calculate_all_web_crippling_conditions(inputs)

    (; t,
    fy,
    panel_depth,
    steel_deck_web_angle_from_horizon,
    steel_deck_outside_radius,
    number_of_webs_per_panel,

    unit_width,

    bearing_width,

    ) = inputs


    #AISI S100-16 Table G5-5

    #One-Flange End

    OFE = calculate_web_crippling_strength(WebCripplingConditionInputs(
        t, fy, panel_depth,
        steel_deck_web_angle_from_horizon,
        steel_deck_outside_radius,
        number_of_webs_per_panel,
        unit_width,
        bearing_width,
        bearing_width / 2,  # N
        4.0,                # C
        0.04,               # C_R
        0.25,               # C_N
        0.025,              # C_h
        1.70,               # Ω_w
        0.90,               # ϕ_w
        0.80,               # ϕ_w_LSD
    ))


    #Two-Flange End

    TFE = calculate_web_crippling_strength(WebCripplingConditionInputs(
        t, fy, panel_depth,
        steel_deck_web_angle_from_horizon,
        steel_deck_outside_radius,
        number_of_webs_per_panel,
        unit_width,
        bearing_width,
        bearing_width / 2,  # N
        9.0,                # C
        0.12,               # C_R
        0.14,               # C_N
        0.040,              # C_h
        1.80,               # Ω_w
        0.85,               # ϕ_w
        0.70,               # ϕ_w_LSD
    ))


    #One-Flange Interior

    OFI = calculate_web_crippling_strength(WebCripplingConditionInputs(
        t, fy, panel_depth,
        steel_deck_web_angle_from_horizon,
        steel_deck_outside_radius,
        number_of_webs_per_panel,
        unit_width,
        bearing_width,
        bearing_width,      # N
        8.0,                # C
        0.10,               # C_R
        0.17,               # C_N
        0.004,              # C_h
        1.75,               # Ω_w
        0.85,               # ϕ_w
        0.75,               # ϕ_w_LSD
    ))


    #Two-Flange Interior

    TFI = calculate_web_crippling_strength(WebCripplingConditionInputs(
        t, fy, panel_depth,
        steel_deck_web_angle_from_horizon,
        steel_deck_outside_radius,
        number_of_webs_per_panel,
        unit_width,
        bearing_width,
        bearing_width,      # N
        10.0,               # C
        0.11,               # C_R
        0.21,               # C_N
        0.02,               # C_h
        1.75,               # Ω_w
        0.85,               # ϕ_w
        0.75,               # ϕ_w_LSD
    ))


    outputs = AllWebCripplingConditions(

            OFE,
            TFE,
            OFI,
            TFI

    )

    return outputs

end
