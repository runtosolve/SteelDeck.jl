

function calculate_all_bare_deck_properties(t, panel_depth, unit_width, fy, E, half_wavelengths, cross_section)


    # #Discretize cross-section:
    # cross_section = CrossSectionGeometry.generate_thin_walled(L, θ, n, r, n_r);


    # #Remove left end leg:
    # cross_section = cross_section[(n[1] + Int(n_r[1]/2) + 1):end];


    # #Scale cross-section by 2 since scale is applied in drawing:
    # cross_section = cross_section .* 2.0;  #hard coded!!!!!

    # #Convert dimensions from mm to inches:
    # cross_section = cross_section ./ 25.4;


    bare_deck_properties_all = Vector{SteelDeck.BareDeckOutputs}(undef, length(t))

    for i in eachindex(t)

            inputs = SteelDeck.BareDeckInputs(
            
            t[i],

            # L,
            # θ,
            # n,
            # r,
            # n_r,
            cross_section, 

            E,
            fy, 

            panel_depth, 
            unit_width, 

            half_wavelengths
            
            ) 


            bare_deck_properties_all[i] = SteelDeck.calculate_bare_deck_properties(inputs)

    end


    return bare_deck_properties_all


end



function calculate_all_bare_shear_properties(t, Fy, E, μ, steel_deck_depth, steel_deck_web_angle_from_horizon, steel_deck_outside_radius, unit_width, number_of_webs_per_panel)

    bare_shear_all = Vector{SteelDeck.BareShearOutputs}(undef, length(t))

    kv = 5.34 

    for i in eachindex(t)

        inputs = SteelDeck.BareShearInputs(
            t[i], 
            Fy, 
            E, 
            μ,
            kv, 

            steel_deck_depth, 
            steel_deck_web_angle_from_horizon,
            steel_deck_outside_radius, 
            number_of_webs_per_panel,

            unit_width 

        )

        bare_shear_all[i] = SteelDeck.calculate_bare_shear_strength(inputs)

    end


    return bare_shear_all

end


function calculate_all_web_crippling_properties(t, Fy, steel_deck_depth, steel_deck_web_angle_from_horizon, steel_deck_outside_radius, bearing_width, unit_width, number_of_webs_per_panel, design_code)

    web_crippling_all = Vector{SteelDeck.AllWebCripplingConditions}(undef, length(t))

    for i in eachindex(t)

        inputs = SteelDeck.WebCripplingInputs(
            t[i],
        steel_deck_depth,
        steel_deck_web_angle_from_horizon,
        steel_deck_outside_radius,
        number_of_webs_per_panel,

        unit_width, 
        

        bearing_width,
        [],

        Fy,

        [], 
        [], 
        [],
        [],

        [],
        [],
        [],

        design_code
        )

        web_crippling_all[i] = SteelDeck.calculate_all_web_crippling_conditions(inputs)

    end

    return web_crippling_all

end




function calculate_all_construction_spans(bare_deck_properties_all, web_crippling_properties_all, bare_shear_properties_all, total_slab_depth, t, deck_flute_widths, dd, panel_width, number_empty_flutes, γc, γs, wlc, Plc, wcdl, E)

    construction_spans_all = OrderedDict{String, OrderedDict}()

    Δ_over_L_limit = 1/180 

    for i in eachindex(t)

        bare_deck_properties = bare_deck_properties_all[i]
        web_crippling_properties = web_crippling_properties_all[i]
        bare_shear_properties = bare_shear_properties_all[i]

        unit_width = 12.0 #in. 

        construction_spans = OrderedDict{String,SteelDeck.ConstructionSpanOutputs}()

        for j in eachindex(total_slab_depth)

            h = total_slab_depth[j]

            inputs = SteelDeck.ConstructionSpanInputs(

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

            )

            construction_spans[string(total_slab_depth[j]) * "inch slab"] = SteelDeck.calculate_construction_spans(inputs)

        end

        construction_spans_all["t = " * string(t[i])] = construction_spans

    end

    return construction_spans_all

end




function calculate_composite_flexural_properties(bare_deck_properties_all, total_slab_depth, Fy, panel_width, Dw, ph, deck_rib_widths, deck_rib_pitch, dd, Es, Wc, fc, β1)

    unit_width = 12.0 #in.

    composite_flexural_all = OrderedDict{String, OrderedDict}()

    for i in eachindex(bare_deck_properties_all)

        composite_flexural = OrderedDict{String,SteelDeck.CompositeFlexuralOutputs}()

        bare_deck_properties = bare_deck_properties_all[i]

        t = bare_deck_properties_all[i].inputs.t

        yc_steel = bare_deck_properties.section_properties.yc
        Ixx_steel_whole_panel = bare_deck_properties.section_properties.Ixx
        As_whole_panel = bare_deck_properties.section_properties.A

        for j in eachindex(total_slab_depth)

            h = total_slab_depth[j] #in. slab depth 

            inputs = SteelDeck.CompositeFlexuralInputs(

                panel_width,
                dd, 
                Dw, 
                ph, 
                As_whole_panel, 
                yc_steel, 
                Ixx_steel_whole_panel, 
                deck_rib_widths,
                deck_rib_pitch,

                h, 

                unit_width, 

                Fy,

                Es, 
                Wc, 
                fc, 
                β1

                ) 


                composite_flexural[string(total_slab_depth[j]) * "inch slab"] = SteelDeck.calculate_composite_flexural_properties(inputs)


        end

        composite_flexural_all["t = " * string(t)] = composite_flexural

    end

    return composite_flexural_all

end



function calculate_all_one_way_shear_strength(t, steel_deck_depth, steel_deck_web_angle_from_horizon, steel_deck_radius, steel_deck_trough_width, number_of_troughs_in_panel, unit_width, total_slab_depth, Fy, E, μ, λ, fc)

    kv = 5.34

    one_way_shear_all = OrderedDict{String, OrderedDict}()

    for i in eachindex(t)

        one_way_shear = OrderedDict{String,SteelDeck.OneWayShearOutputs}()

        for j in eachindex(total_slab_depth)

            inputs = SteelDeck.OneWayShearInputs(

                t[i], 
                steel_deck_depth,
                steel_deck_web_angle_from_horizon,
                steel_deck_radius,
                steel_deck_trough_width, 
                number_of_troughs_in_panel,
                unit_width,

                total_slab_depth[j], 

                Fy,
                E, 
                μ,

                kv,

                λ,
                fc

            )

            one_way_shear[string(total_slab_depth[j]) * "inch slab"] = SteelDeck.calculate_one_way_composite_shear_strength(inputs)

        end

        one_way_shear_all["t = " * string(t[i])] = one_way_shear

    end

    return one_way_shear_all

end





function calculate_all_superimposed_dead_load(composite_deck_flexural_properties_all, construction_spans_all, one_way_shear_all, total_slab_depth, t, span_lengths, E)

    Mn_ASD_all = []
    Vn_ASD_all = []
    W1_concrete_all = []
    W1_steel_all = []
    Id_all = []

    Δ_over_L_limit = 1 / 360 

    for i in eachindex(composite_deck_flexural_properties_all)

        for j in keys(composite_deck_flexural_properties_all[i])

            if isnothing(composite_deck_flexural_properties_all[i][j].Mno_ASD)

                push!(Mn_ASD_all, composite_deck_flexural_properties_all[i][j].Mro_ASD)

            else

                push!(Mn_ASD_all, composite_deck_flexural_properties_all[i][j].Mno_ASD)

            end

            push!(W1_concrete_all, construction_spans_all[i][j].W1_concrete)
            push!(W1_steel_all, construction_spans_all[i][j].W1_steel)
            push!(Vn_ASD_all, one_way_shear_all[i][j].aVn_unit)
            push!(Id_all, composite_deck_flexural_properties_all[i][j].Id)

        end


    end

    num_panel_thicknesses = length(composite_deck_flexural_properties_all);
    num_slab_thicknesses = length(collect(keys(composite_deck_flexural_properties_all[collect(keys(composite_deck_flexural_properties_all))[1]])));

    #+ echo=false
    Mn_ASD_all = reshape(Mn_ASD_all, (num_slab_thicknesses, num_panel_thicknesses));
    Vn_ASD_all = reshape(Vn_ASD_all, (num_slab_thicknesses, num_panel_thicknesses));
    W1_concrete_all = reshape(W1_concrete_all, (num_slab_thicknesses, num_panel_thicknesses));
    W1_steel_all = reshape(W1_steel_all, (num_slab_thicknesses, num_panel_thicknesses));
    Id_all = reshape(Id_all, (num_slab_thicknesses, num_panel_thicknesses));



    all_superimposed_load_results = Array{SteelDeck.SuperimposedLoadOutputs}(undef, num_slab_thicknesses, length(span_lengths), num_panel_thicknesses)

    for i=1:num_slab_thicknesses

        for j= 1:num_panel_thicknesses

            Mn_ASD = Mn_ASD_all[i, j]
            Vn_ASD = Vn_ASD_all[i, j]
            W1_steel = W1_steel_all[i, j]
            W1_concrete = W1_concrete_all[i, j]
            Id = Id_all[i, j]

            panel_thickness = t[j]

            for k in eachindex(span_lengths)

            
                span_length = span_lengths[k]
            
                inputs = SteelDeck.SuperimposedLoadInputs(
                    E,
                    Id, 
                    
                    Mn_ASD,
                    Vn_ASD, 

                    panel_thickness,
                    total_slab_depth[i],
                    span_length, 
                    
                    W1_steel, 
                    W1_concrete, 
                    
                    Δ_over_L_limit)


                all_superimposed_load_results[i, k, j] = SteelDeck.calculate_superimposed_load(inputs)

            end

        end

    end

    return all_superimposed_load_results

end


function generate_bare_deck_section_properties_table(t_panel, t_gage, bare_deck_all)

    Ip = [bare_deck_all[i].I_eff_pos_unit for i in eachindex(bare_deck_all)];
    In = [bare_deck_all[i].I_eff_neg_unit for i in eachindex(bare_deck_all)];
    Sp = [bare_deck_all[i].S_pos_eff_unit for i in eachindex(bare_deck_all)];
    Sn = [bare_deck_all[i].S_neg_eff_unit for i in eachindex(bare_deck_all)];

    table = DataFrame();

    table.gage = t_gage;
    table.design_thickness = t_panel;
    table.Ip = round.(Ip, digits=4);
    table.In = round.(In, digits=4); 
    table.Sp = round.(Sp, digits =4);
    table.Sn = round.(Sn, digits=4);

    return table 

end


function generate_bare_deck_web_crippling_table(web_crippling_all, t_gage)

    aPn_OFE = [web_crippling_all[i].OFE.aPn_unit for i in eachindex(web_crippling_all)];
    aPn_TFE = [web_crippling_all[i].TFE.aPn_unit for i in eachindex(web_crippling_all)];

    aPn_OFI = [web_crippling_all[i].OFI.aPn_unit for i in eachindex(web_crippling_all)];
    aPn_TFI = [web_crippling_all[i].TFI.aPn_unit for i in eachindex(web_crippling_all)];

    table = DataFrame();
    table.gage = t_gage;
    table.OFE = Int.(round.(aPn_OFE, digits=0));
    table.OFI = Int.(round.(aPn_OFI, digits=0));
    table.TFE = Int.(round.(aPn_TFE, digits=0));
    table.TFI = Int.(round.(aPn_TFI, digits=0));

    return table 

end


function generate_bare_deck_allowable_shear_strength_table(bare_shear_all, t_gage)

    aVn_ASD = [bare_shear_all[i].aVn_unit_ASD for i in eachindex(bare_shear_all)];
    # aVn_LRFD = [bare_shear_all[i].aVn_unit_LRFD for i in eachindex(bare_shear_all)];

    table = DataFrame();

    table.gage = t_gage;
    table.aVn_ASD = Int.(round.(aVn_ASD, digits=0));


    return table 

end


function generate_bare_deck_allowable_flexural_strength_table(bare_deck_all, t_gage)

    aMn_P_ASD = [bare_deck_all[i].aMnℓ_pos_unit_ASD for i in eachindex(bare_deck_all)];
    aMn_N_ASD = [bare_deck_all[i].aMnℓ_neg_unit_ASD for i in eachindex(bare_deck_all)];

    table = DataFrame();

    table.gage = t_gage;
    table.aMn_P_ASD = Int.(round.([aMn_P_ASD[i] for i in eachindex(aMn_P_ASD)]));
    table.aMn_N_ASD = Int.(round.([aMn_N_ASD[i] for i in eachindex(aMn_N_ASD)]));

    return table 

end


function generate_composite_deck_flexural_property_tables(composite_flexural_all)

    Iu = [];
    Icr = [];
    Id = [];
    slab_depth = [];
    My = [];
    aMn = [];

    for i in keys(composite_flexural_all)

         for j in keys(composite_flexural_all[i])

              push!(Iu, composite_flexural_all[i][j].Iu)
              push!(Icr, composite_flexural_all[i][j].Icr)
              push!(Id, composite_flexural_all[i][j].Id)
              push!(slab_depth, composite_flexural_all[i][j].inputs.h)
              push!(My, composite_flexural_all[i][j].My)

              if isnothing(composite_flexural_all[i][j].Mno_ASD)

                   push!(aMn, composite_flexural_all[i][j].Mro_ASD)

              else

                   push!(aMn, composite_flexural_all[i][j].Mno_ASD)

              end

         end

    end;

    num_panel_thicknesses = length(composite_flexural_all);
    num_slab_thicknesses = length(collect(keys(composite_flexural_all[collect(keys(composite_flexural_all))[1]])));


    Iu = reshape(Iu, (num_slab_thicknesses, num_panel_thicknesses));
    Icr = reshape(Icr, (num_slab_thicknesses, num_panel_thicknesses));
    Id = reshape(Id, (num_slab_thicknesses, num_panel_thicknesses));
    slab_depth = reshape(slab_depth, (num_slab_thicknesses, num_panel_thicknesses));
    My = reshape(My, (num_slab_thicknesses, num_panel_thicknesses));
    aMn = reshape(aMn, (num_slab_thicknesses, num_panel_thicknesses));

    table = DataFrame();
    table.slab_depth = slab_depth[:, 1];
    table.Iu_gage_18 = round.(Iu[:, 2], digits=2);
    table.Iu_gage_20 = round.(Iu[:, 3], digits=2);
    table.Iu_gage_22 = round.(Iu[:, 4], digits=2);


    table_uncracked = table 

    table = DataFrame();
    table.slab_depth = slab_depth[:, 1];
    table.Icr_gage_18 = round.(Icr[:, 2], digits=2);
    table.Icr_gage_20 = round.(Icr[:, 3], digits=2);
    table.Icr_gage_22 = round.(Icr[:, 4], digits=2);

    table_cracked = table

    table = DataFrame();
    table.slab_depth = slab_depth[:, 1];
    table.Id_gage_18 = round.(Id[:, 2], digits=2);
    table.Id_gage_20 = round.(Id[:, 3], digits=2);
    table.Id_gage_22 = round.(Id[:, 4], digits=2);

    table_service = table 

    table = DataFrame();
    table.slab_depth = slab_depth[:, 1];
    table.aMn_gage_18 = Int.(round.(aMn[:, 2]));
    table.aMn_gage_20 = Int.(round.(aMn[:, 3]));
    table.aMn_gage_22 = Int.(round.(aMn[:, 4]));

    table_allowable = table 

    return table_uncracked, table_cracked, table_service, table_allowable

end


function write_feet_inches(ℓ)

    remainder = floor(Int, (ℓ - floor(Int, ℓ)) * 12) 
    measurement = string(floor(Int, ℓ)) * "'-" * string(remainder) * "\""

    return measurement

end;


function collect_spans(ℓ_min, slab_depth)

    simple_span = []
    double_span = []
    triple_span = []

    for i=1:size(slab_depth)[1]
         simple_span = [simple_span; write_feet_inches.([ℓ_min[i, j]["simple span"] for j=1:4] ./ 12)];
         double_span = [double_span; write_feet_inches.([ℓ_min[i, j]["double span"] for j=1:4] ./ 12)];
         triple_span = [triple_span; write_feet_inches.([ℓ_min[i, j]["triple span"] for j=1:4] ./ 12)];
    end

    return simple_span, double_span, triple_span

end;



function generate_construction_span_table(construction_spans_all, t_gage)

  
    ℓ_min = [];
    slab_depth = [];

    for i in keys(construction_spans_all)

         for j in keys(construction_spans_all[i])

              push!(ℓ_min, construction_spans_all[i][j].ℓ_min)
              push!(slab_depth, construction_spans_all[i][j].inputs.h)

         end

    end;

    num_panel_thicknesses = length(construction_spans_all);
    num_slab_thicknesses = length(collect(keys(construction_spans_all[collect(keys(construction_spans_all))[1]])));

   
    ℓ_min = reshape(ℓ_min, (num_slab_thicknesses, num_panel_thicknesses));
    slab_depth = reshape(slab_depth, (num_slab_thicknesses, num_panel_thicknesses));
    gages = repeat(t_gage, num_slab_thicknesses);


 
    simple_span, double_span, triple_span = collect_spans(ℓ_min, slab_depth);

    slab_depth = vec(slab_depth[:, 1:4]');;


    table = DataFrame();
    table.slab_depth = slab_depth;
    table.gages = gages;
    table.ℓ_simple_span = simple_span;
    table.ℓ_double_span = double_span;
    table.ℓ_triple_span = triple_span;

    return table 

end





function generate_composite_deck_one_way_shear_strength_table(one_way_shear_all)

    #+ echo=false
    aVn = [];
    slab_depth = [];

    for i in keys(one_way_shear_all)

         for j in keys(one_way_shear_all[i])

              push!(aVn, one_way_shear_all[i][j].aVn_unit)
              push!(slab_depth, one_way_shear_all[i][j].inputs.total_slab_depth)

         end

    end;
 
    num_panel_thicknesses = length(one_way_shear_all);
    num_slab_thicknesses = length(collect(keys(one_way_shear_all[collect(keys(one_way_shear_all))[1]])));

    aVn = reshape(aVn, (num_slab_thicknesses, num_panel_thicknesses));
    slab_depth = reshape(slab_depth, (num_slab_thicknesses, num_panel_thicknesses));

 
    table = DataFrame();
    table.slab_depth = slab_depth[:, 1];
    table.aVn_gage_18 = Int.(round.(aVn[:, 2]));
    table.aVn_gage_20 = Int.(round.(aVn[:, 3]));
    table.aVn_gage_22 = Int.(round.(aVn[:, 4]));

    return table 


end



function generate_superimposed_dead_load_table(superimposed_loads_all, total_slab_depth, slab_spans, gage_index)

    superimposed_load = [];
    deck_slab_weight = [];
    deck_slab_depth = [];


    # i = 2;
    i = gage_index
    for j in eachindex(total_slab_depth)
         for k in eachindex(slab_spans)

              push!(superimposed_load, superimposed_loads_all[j, k, i].superimposed_load)
              push!(deck_slab_weight, superimposed_loads_all[j, k, i].deck_slab_weight)
              push!(deck_slab_depth, superimposed_loads_all[j, k, i].inputs.total_slab_depth)

         end

    end


    superimposed_load = reshape(superimposed_load, (length(slab_spans), length(total_slab_depth)))';
    deck_slab_weight = reshape(deck_slab_weight, (length(slab_spans), length(total_slab_depth)))';
    deck_slab_depth = reshape(deck_slab_depth, (length(slab_spans), length(total_slab_depth)))';


    span_labels = [write_feet_inches(slab_spans[i]) for i in eachindex(slab_spans)];


    table = DataFrame([deck_slab_depth[:, 1] Int.(round.([deck_slab_weight[:, 1] superimposed_load]))], ["deck_slab_depth"; "deck_slab_weight"; span_labels]);
    header = (["slab depth"; "slab weight"; span_labels], ["in"; fill("psf", 15)]);

    return table 

end
