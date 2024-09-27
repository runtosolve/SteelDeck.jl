module SteelDeck

using CUFSM, SectionProperties, AISIS100


struct BareGeometry

    # L::Vector{Float64}
    # θ::Vector{Float64}
    # n::Vector{Int}
    # r::Vector{Float64}
    # n_r::Vector{Int}
    t::Float64

    cross_section::Vector{Vector{Float64}}
    X::Vector{Float64}
    Y::Vector{Float64}

end


struct BareProperties

    section_properties::CUFSM.SectionPropertiesObject

    fy::Float64
    E::Float64 

    panel_depth::Float64
    unit_width::Float64

    local_buckling_pos::CUFSM.Model
    local_buckling_neg::CUFSM.Model
    
    Mcrℓ_pos::Float64
    Mcrℓ_neg::Float64 

    Mcrℓ_pos_unit::Float64
    Mcrℓ_neg_unit::Float64 

    y_top::Float64
    y_bottom::Float64

    Ixx::Float64
    Ixx_unit::Float64

    S_pos_unit::Float64
    S_neg_unit::Float64

    My_pos_unit::Float64
    My_neg_unit::Float64
    My_unit::Float64

    Mnℓ_pos_unit::Float64
    Mnℓ_neg_unit::Float64

    aMnℓ_pos_unit_ASD::Float64
    aMnℓ_neg_unit_ASD::Float64

    aMnℓ_pos_unit_LRFD::Float64
    aMnℓ_neg_unit_LRFD::Float64
    
    Md_pos_unit::Float64
    Md_neg_unit::Float64

    I_eff_pos_unit::Float64
    I_eff_neg_unit::Float64

    S_pos_eff_unit::Float64
    S_neg_eff_unit::Float64
 
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

function calculate_bare_properties(inputs)

    #Unpack inputs:
    cross_section, t, fy, half_wavelengths, E, panel_depth, unit_width = inputs

    #Find number of cross-section elements:
    num_elem = size(cross_section)[1] - 1;

    #Calculate section properties:
    center = cross_section
    section_properties = SectionProperties.open_thin_walled(center, t * ones(Float64, num_elem))


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

    y_top = panel_depth - section_properties.yc 
    y_bottom = section_properties.yc 

    S_pos_unit = Ixx_unit / y_top
    S_neg_unit = Ixx_unit / y_bottom

    #Calculate yield moment per unit width:
    My_pos_unit = fy * S_pos_unit
    My_neg_unit = fy * S_neg_unit
    My_unit = minimum([My_pos_unit, My_neg_unit])


    #Calculate ASD Mnl:
    design_code = "AISI S100-16 ASD"
    Mnℓ_pos_unit, aMnℓ_pos_unit_ASD = AISIS100.v16S3.f321(My_unit, Mcrℓ_pos_unit, design_code)
    Mnℓ_neg_unit, aMnℓ_neg_unit_ASD = AISIS100.v16S3.f321(My_unit, Mcrℓ_neg_unit, design_code)

    #Calculate LRFD Mnl:
    design_code = "AISI S100-16 LRFD"
    Mnℓ_pos_unit, aMnℓ_pos_unit_LRFD = AISIS100.v16S3.f321(My_unit, Mcrℓ_pos_unit, design_code)
    Mnℓ_neg_unit, aMnℓ_neg_unit_LRFD = AISIS100.v16S3.f321(My_unit, Mcrℓ_neg_unit, design_code)


    #Calculate the effective moment of inertia and section moduli:
    Md_pos_unit, not_used = AISIS100.v16S3.f321(aMnℓ_pos_unit_ASD, Mcrℓ_pos_unit, design_code)
    Ixx_eff_pos_unit = AISIS100.v16S3.l21(Md_pos_unit, aMnℓ_pos_unit_ASD, Ixx_unit)
    S_pos_eff_unit = Ixx_eff_pos_unit / y_top

    Md_neg_unit, not_used = AISIS100.v16S3.f321(aMnℓ_neg_unit_ASD, Mcrℓ_neg_unit, design_code)
    Ixx_eff_neg_unit = AISIS100.v16S3.l21(Md_neg_unit, aMnℓ_neg_unit_ASD, Ixx_unit)
    S_neg_eff_unit = Ixx_eff_neg_unit / y_bottom

    bare_properties = BareProperties(
                    section_properties,
                    fy,
                    E,

                    panel_depth,
                    unit_width,

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

                    My_pos_unit,
                    My_neg_unit,
                    My_unit,

                    Mnℓ_pos_unit,
                    Mnℓ_neg_unit,

                    aMnℓ_pos_unit_ASD,
                    aMnℓ_neg_unit_ASD,

                    aMnℓ_pos_unit_LRFD,
                    aMnℓ_neg_unit_LRFD,

                    Md_pos_unit,
                    Md_neg_unit,

                    Ixx_eff_pos_unit,
                    Ixx_eff_neg_unit,

                    S_pos_eff_unit,
                    S_neg_eff_unit)

    return bare_properties 

end



function calculate_all_bare_properties(inputs)

    #t should be a vector here 
    cross_section, t, fy, half_wavelengths, E, panel_depth, unit_width = inputs

    bare_geometry_all = Vector{SteelDeck.BareGeometry}(undef, length(t))
    bare_properties_all = Vector{SteelDeck.BareProperties}(undef, length(t))

    for i in eachindex(t)

        #Zero cross-section considering thickness:

        X = [cross_section[i][1] for i in eachindex(cross_section)];
        Y = [cross_section[i][2] for i in eachindex(cross_section)];

        ΔX = -minimum(X)
        ΔY = -minimum(Y) + t[i] / 2

        X .+= ΔX
        Y .+= ΔY

        cross_section = [[X[i], Y[i]] for i in eachindex(cross_section)]

        bare_geometry_all[i] = SteelDeck.BareGeometry(

            # L,
            # θ,
            # n,
            # r,
            # n_r,
            t[i],

            cross_section,
            X,
            Y
        )

        inputs = (cross_section, t[i], fy, half_wavelengths, E, panel_depth, unit_width)

        cross_section, t[i], fy, half_wavelengths, E, panel_depth, unit_width = inputs

        bare_properties_all[i] = SteelDeck.calculate_bare_properties(inputs)

    end

    return bare_geometry_all, bare_properties_all

end

end # module SteelDeck
