


struct CompositeFlexuralInputs

    panel_width
    dd 
    Dw 
    ph 
    As_whole_panel 
    yc_steel 
    Ixx_steel_whole_panel 
    deck_rib_widths
    deck_rib_pitch

    h 

    unit_width 

    Fy

    Es
    Wc
    fc
    β1
    embossement_type
    t
    ps1
    ps2
    design_method

    L

end

function CompositeFlexuralInputs(panel_width, dd, Dw, ph, As_whole_panel, yc_steel, Ixx_steel_whole_panel, deck_rib_widths, deck_rib_pitch, h, unit_width, Fy, Es, Wc, fc, β1, embossement_type, t, design_method, L; ps1=0, ps2=0)
    CompositeFlexuralInputs(panel_width, dd, Dw, ph, As_whole_panel, yc_steel, Ixx_steel_whole_panel, deck_rib_widths, deck_rib_pitch, h, unit_width, Fy, Es, Wc, fc, β1, embossement_type, t, ps1, ps2, design_method, L)
end



struct CompositeFlexuralOutputs 

    inputs 

    b 
    Ec 
    n 
    hc 
    d 
    Wr
    Cs

    ycc_uncracked
    ycs_uncracked
    Iu
    Id

    ρ 

    ycc 
    ycs 
    
    As 
    Isf 
    
    Icr 

    K1 
    K3 
    K 

    My 

    Mn
    c_over_d
    c_over_d_b

   
end


function calculate_composite_flexural_properties(inputs)

    (;

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
    β1,
    embossement_type,
    t,
    ps1,
    ps2,
    design_method,

    L,

    ) = inputs


    Wr = mean(deck_rib_widths)

    Cs = deck_rib_pitch

    Isf = Ixx_steel_whole_panel * unit_width / panel_width

    b = unit_width 

    Ec = Wc^1.5 * (fc/1000)^0.5 * 1000

    n = Es / Ec

    hc = h - dd 

    d = h - yc_steel

    As = As_whole_panel * unit_width / panel_width

    ρ = As / (b * d)


    ### uncracked 
    
    
    ycc_uncracked = SDIComposite.C2022.EqC_F2_3_3(b, hc, n, As, d, Wr, dd, h, Cs)

    ycs_uncracked = d - ycc_uncracked 

    Iu = SDIComposite.C2022.EqC_F2_3_4(b, h, hc, n, ycc_uncracked, Isf, As, ycs_uncracked, Wr, dd, Cs)

  

    #### cracked 

    ycc = SDIComposite.C2022.EqC_F2_3_1(d, ρ, n, hc)

    ycs = d - ycc 

    Icr = SDIComposite.C2022.EqC_F2_3_2(b, n, ycc, As, ycs, Isf)

    Id = SDIComposite.C2022.EqC_F2_3_5(Iu, Icr)

    if embossement_type == "type_1"
            K1 = SDIComposite.C2022.EqF3_2__7_type_1(Dw, ph)
    elseif embossement_type == "type_2"
            K1 = SDIComposite.C2022.EqF3_2__8_type_2(t, Dw, ph)
    elseif embossement_type == "type_3"
            K11 = SDIComposite.C2022.EqF3_2__7_type_1(Dw, ph)
            K12 = SDIComposite.C2022.EqF3_2__8_type_2(t, Dw, ph)
            K1 = SDIComposite.C2022.EqF3_2__9_type_3(K11, K12, ps1, ps2)
    else
        error("Invalid embossement type. Must be type_1, type_2, or type_3.")
    end

    K3 = 1.4


    c_over_d = SDIComposite.C2022.EqC_F2_2__1(As, Fy, fc, d, b, β1)

    c_over_d_b = SDIComposite.C2022.EqC_F2_2__2(h, dd, Fy, Es, d)

    My = SDIComposite.C2022.EqF3_2__11(Fy, Icr, h, ycc)

    K = SDIComposite.C2022.EqF_3_2__6(K1, K3)


    Mn = SDIComposite.C2022.EqF3_2__10(K, My, L, design_method)

    outputs = CompositeFlexuralOutputs(

        inputs,

        b,
        Ec,
        n,
        hc,
        d,
        Wr,
        Cs,

        ycc_uncracked,
        ycs_uncracked,
        Iu,
        Id,

        ρ,

        ycc,
        ycs,

        As,
        Isf,

        Icr,

        K1,
        K3,
        K,

        My,

        Mn,
        c_over_d,
        c_over_d_b,

    )

     return outputs 

end