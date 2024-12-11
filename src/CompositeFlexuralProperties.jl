


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

    c_over_d 
    c_over_d_b

    Mno_ASD
    Mno_LRFD

    ϵcu
    m 
    c
    Mro_ASD 
    Mro_LRFD 
   
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
    β1

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
    
    
    ycc_uncracked = SDIComposite.C2017.EqA5__3(b, hc, n, As, d, Wr, dd, h, Cs)

    ycs_uncracked = d - ycc_uncracked 

    Iu = SDIComposite.C2017.EqA5__4(b, h, hc, n, ycc_uncracked, Isf, As, ycs_uncracked, Wr, dd, Cs)

  

    #### cracked 

    ycc = SDIComposite.C2017.EqA5__1(d, ρ, n, hc)

    ycs = d - ycc 

    Icr = SDIComposite.C2017.EqA5__2(b, n, ycc, As, ycs, Isf)

    Id = SDIComposite.C2017.EqA5__5(Iu, Icr)

    K1 = SDIComposite.C2017.EqA2__12(Dw, ph)
    K3 = 1.4


    c_over_d = SDIComposite.C2017.EqA2__5(As, Fy, fc, d, b, β1)

        c_over_d_b = SDIComposite.C2017.EqA2__6(h, dd, Fy, Es, d)

        My = SDIComposite.C2017.EqA2__9(Fy, Icr, h, ycc)

        K = SDIComposite.C2017.EqA2__10(K1, K3)

    if c_over_d < c_over_d_b
        Mno_ASD =  SDIComposite.C2017.EqA2__15(K, My)
        Mno_LRFD =  SDIComposite.C2017.EqA2__16(K, My)

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
        
            c_over_d, 
            c_over_d_b,
        
            Mno_ASD,
            Mno_LRFD,
        
            nothing,
            nothing, 
            nothing,
            nothing, 
            nothing
        ) 
        



    else 

        ϵcu = 0.003

        m = SDIComposite.C2017.EqA2__20(Es, ϵcu, fc, β1)

        c = SDIComposite.C2017.EqA2__18(d, ρ, m)

        Mro_LRFD = SDIComposite.C2017.EqA2__17a(fc, b, β1, c, d, My, K)

        Mro_ASD = SDIComposite.C2017.EqA2__17b(fc, b, β1, c, d, My, K)


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
        
            c_over_d, 
            c_over_d_b,
        
            nothing,
            nothing,
        
            ϵcu,
            m, 
            c,
            Mro_ASD, 
            Mro_LRFD
        ) 

    end

    return outputs 

end