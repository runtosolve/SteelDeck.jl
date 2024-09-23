module SteelDeck

using CUFSM 

function calculate_bare_deck_local_buckling(cross_section, t_panel, E, lengths)

    local_buckling_pos = Vector{CUFSM.Model}(undef, length(t_panel))
    Mcrl_pos = Vector{Float64}(undef, length(t_panel))

    for i in eachindex(t_panel)

        cross_section = center_all[i]
        x_center = [cross_section[i][1] for i in eachindex(cross_section)];
        y_center = [cross_section[i][2] for i in eachindex(cross_section)];

        t = t_panel[i] * ones(Float64, num_elem) 

        lengths = 4.0:1.0:10.0
        E = 29500.0
        ν = 0.30
        P = 0.0
        Mxx = 1.0
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        springs = []
        constraints = []
        supports = [[1 0 0 0 0], [num_elem+1, 0, 0, 0, 0]]
        neigs = 1

        model = CUFSM.Tools.open_section_analysis(x_center, y_center, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, supports, neigs)

        local_buckling_pos[i] = model 
        eig = 1
        Mcrl_pos[i] = minimum(CUFSM.Tools.get_load_factor(model, eig))

    end

    return Mcrl_pos 

end



end # module SteelDeck
