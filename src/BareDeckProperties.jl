struct BareGeometry

    t::Float64

    cross_section::Vector{Vector{Float64}}
    X::Vector{Float64}
    Y::Vector{Float64}

end


struct BareDeckInputs

    t

    # L
    # θ
    # n
    # r
    # n_r
    cross_section

    E

    panel_depth
    unit_width

    half_wavelengths

    design_method

end


Base.@kwdef struct BareDeckOutputs

    inputs = nothing

    section_properties = nothing

    local_buckling_pos = nothing
    local_buckling_neg = nothing

    Mcrℓ_pos = nothing
    Mcrℓ_neg = nothing

    Mcrℓ_pos_unit = nothing
    Mcrℓ_neg_unit = nothing

    y_top = nothing
    y_bottom = nothing

    Ixx = nothing
    Ixx_unit = nothing

    S_pos_unit = nothing
    S_neg_unit = nothing
    Sp_unit = nothing


    t = nothing

    fy = nothing

    My_pos_unit = nothing
    My_neg_unit = nothing
    My_unit = nothing
    M_plastic_unit = nothing

    Mnℓ_pos_unit = nothing
    Mnℓ_neg_unit = nothing

    aMnℓ_pos_unit = nothing
    aMnℓ_neg_unit = nothing

    Md_pos_unit = nothing
    Md_neg_unit = nothing

    I_eff_pos_unit = nothing
    I_eff_neg_unit = nothing

    S_pos_eff_unit = nothing
    S_neg_eff_unit = nothing

end






function calculate_bare_deck_local_buckling(cross_section, t, E, lengths, Mxx)

    num_elem = size(cross_section)[1] - 1

    x_center = [cross_section[i][1] for i in eachindex(cross_section)];
    y_center = [cross_section[i][2] for i in eachindex(cross_section)];

    t = t * ones(Float64, num_elem)

    ν = 0.30
    P = 0.0

    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0

    springs = []
    constraints = []
    supports = [[1 0 0 0 0], [num_elem+1, 0, 0, 0, 0]]
    neigs = 1

    model = CUFSM.Tools.open_section_analysis(x_center, y_center, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, supports, neigs)

    return model

end


function calculate_bare_deck_base_properties(inputs)

    #Unpack inputs:
    (;

    t,

    # L,
    # θ,
    # n,
    # r,
    # n_r,
    cross_section,

    E,

    panel_depth,
    unit_width,

    half_wavelengths

    ) = inputs

    #Find number of cross-section elements:
    num_elem = size(cross_section)[1] - 1;

    #Calculate section properties:
    center = cross_section
    section_properties = SectionProperties.open_thin_walled(center, t * ones(Float64, num_elem))

    #Calculate plastic section properties:
    X = [cross_section[i][1] for i in eachindex(cross_section)]
    Y = [cross_section[i][2] for i in eachindex(cross_section)]

    num_nodes        = length(X)
    node_geometry    = Float64[X Y]
    node_start       = Float64.(collect(1:num_nodes-1))
    node_end         = Float64.(collect(2:num_nodes))
    elem_thicknesses = fill(t, num_nodes-1)
    element_definitions = hcat(node_start, node_end, elem_thicknesses)
    plastic_properties = SectionProperties.calculate_plastic_section_properties(node_geometry, element_definitions, "x")
    Sp = plastic_properties.Z
    Sp_unit = Sp / unit_width

    #Calculate positive local buckling moment:
    lengths = half_wavelengths
    Mxx = +1.0
    model = calculate_bare_deck_local_buckling(cross_section, t, E, lengths, Mxx)

    local_buckling_pos = model

    eig = 1
    Mcrℓ_pos = minimum(CUFSM.Tools.get_load_factor(model, eig))

    #Calculate negative local buckling moment:
    lengths = half_wavelengths
    Mxx = -1.0
    model = calculate_bare_deck_local_buckling(cross_section, t, E, lengths, Mxx)

    local_buckling_neg = model

    eig = 1
    Mcrℓ_neg = minimum(CUFSM.Tools.get_load_factor(model, eig))

    #Calculate Mcrℓ per unit width:
    Mcrℓ_pos_unit = Mcrℓ_pos / unit_width
    Mcrℓ_neg_unit = Mcrℓ_neg / unit_width

    #Calculate gross section moduli per unit width:
    Ixx = section_properties.Ixx
    Ixx_unit = Ixx / unit_width

    y_bottom = section_properties.yc
    y_top = panel_depth - y_bottom

    S_pos_unit = Ixx_unit / y_top
    S_neg_unit = Ixx_unit / y_bottom

    outputs = BareDeckOutputs(;

        inputs,

        section_properties,

        local_buckling_pos,
        local_buckling_neg,

        Mcrℓ_pos,
        Mcrℓ_neg,

        Mcrℓ_pos_unit,
        Mcrℓ_neg_unit,

        y_top,
        y_bottom,

        Ixx,
        Ixx_unit,

        S_pos_unit,
        S_neg_unit,
        Sp_unit,

    )

    return outputs

end


function calculate_bare_deck_properties(base_outputs, fy, S100_version)

    #Unpack base outputs:
    (;

    t,
    Mcrℓ_pos_unit,
    Mcrℓ_neg_unit,
    y_top,
    y_bottom,
    Ixx_unit,
    S_pos_unit,
    S_neg_unit,
    Sp_unit,
    inputs,
    ) = base_outputs

    panel_depth   = inputs.panel_depth
    design_method = inputs.design_method

    #Calculate plastic moment per unit width:
    M_plastic_unit = Sp_unit * fy

    #Calculate yield moment per unit width:
    My_pos_unit = fy * S_pos_unit
    My_neg_unit = fy * S_neg_unit
    My_unit = minimum([My_pos_unit, My_neg_unit])

    ks  = M_plastic_unit / My_unit
    αs  = 1.0
    d   = panel_depth
    My3 = M_plastic_unit - (M_plastic_unit - My_unit) / 9

    if S100_version == "v24"

        #Calculate Mnl (positive):
        βs_pos = max(2*y_top/d, 0.4)
        Mnℓ_pos_unit, aMnℓ_pos_unit = v2024.f32(My_unit, Mcrℓ_pos_unit, ks, αs, βs_pos, My3, design_method)

        #Calculate Mnl (negative):
        βs_neg = max(2*y_bottom/d, 0.4)
        Mnℓ_neg_unit, aMnℓ_neg_unit = v2024.f32(My_unit, Mcrℓ_neg_unit, ks, αs, βs_neg, My3, design_method)

        #Calculate the effective moment of inertia and section moduli:
        Md_pos_unit, _ = v2024.f32(aMnℓ_pos_unit, Mcrℓ_pos_unit, ks, αs, βs_pos, My3, design_method)
        I_eff_pos_unit = AISIS100.v16S3.l21(Md_pos_unit, aMnℓ_pos_unit, Ixx_unit)
        S_pos_eff_unit = I_eff_pos_unit / y_top

        Md_neg_unit, _ = v2024.f32(aMnℓ_neg_unit, Mcrℓ_neg_unit, ks, αs, βs_neg, My3, design_method)
        I_eff_neg_unit = AISIS100.v16S3.l21(Md_neg_unit, aMnℓ_neg_unit, Ixx_unit)
        S_neg_eff_unit = I_eff_neg_unit / y_bottom

    else  # v16

        method_str = "AISI S100-16 $design_method"

        Mnℓ_pos_unit, aMnℓ_pos_unit = AISIS100.v16S3.f321(My_unit, Mcrℓ_pos_unit, method_str)
        Mnℓ_neg_unit, aMnℓ_neg_unit = AISIS100.v16S3.f321(My_unit, Mcrℓ_neg_unit, method_str)

        #Calculate the effective moment of inertia and section moduli:
        Md_pos_unit, _ = AISIS100.v16S3.f321(aMnℓ_pos_unit, Mcrℓ_pos_unit, method_str)
        I_eff_pos_unit = AISIS100.v16S3.l21(Md_pos_unit, aMnℓ_pos_unit, Ixx_unit)
        S_pos_eff_unit = I_eff_pos_unit / y_top

        Md_neg_unit, _ = AISIS100.v16S3.f321(aMnℓ_neg_unit, Mcrℓ_neg_unit, method_str)
        I_eff_neg_unit = AISIS100.v16S3.l21(Md_neg_unit, aMnℓ_neg_unit, Ixx_unit)
        S_neg_eff_unit = I_eff_neg_unit / y_bottom

    end

    outputs = BareDeckOutputs(;

        inputs            = base_outputs.inputs,
        section_properties = base_outputs.section_properties,
        local_buckling_pos = base_outputs.local_buckling_pos,
        local_buckling_neg = base_outputs.local_buckling_neg,
        Mcrℓ_pos          = base_outputs.Mcrℓ_pos,
        Mcrℓ_neg          = base_outputs.Mcrℓ_neg,
        Mcrℓ_pos_unit     = base_outputs.Mcrℓ_pos_unit,
        Mcrℓ_neg_unit     = base_outputs.Mcrℓ_neg_unit,
        y_top             = base_outputs.y_top,
        y_bottom          = base_outputs.y_bottom,
        Ixx               = base_outputs.Ixx,
        Ixx_unit          = base_outputs.Ixx_unit,
        S_pos_unit        = base_outputs.S_pos_unit,
        S_neg_unit        = base_outputs.S_neg_unit,
        Sp_unit           = base_outputs.Sp_unit,


        t,

        fy,

        My_pos_unit,
        My_neg_unit,
        My_unit,
        M_plastic_unit,

        Mnℓ_pos_unit,
        Mnℓ_neg_unit,

        aMnℓ_pos_unit,
        aMnℓ_neg_unit,

        Md_pos_unit,
        Md_neg_unit,

        I_eff_pos_unit,
        I_eff_neg_unit,

        S_pos_eff_unit,
        S_neg_eff_unit,

    )

    return outputs

end