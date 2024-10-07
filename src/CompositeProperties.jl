



function calculate_composite_flexural_strength(Es, Fy, fc, β1, d, dd, b, h, ρ, As, Icr, ycc)


    c_over_d = SDIComposite.C2017.EqA3__3(As, Fy, fc, d, b, β1)

    c_over_d_b = SDIComposite.C2017.EqA3__4(h, dd, Fy, Es, d)

    My = (Fy * Icr) / (h - ycc)

    if c_over_d < c_over_d_b

        Mr_LRFD =  SDIComposite.C2017.Eq_A3__5a(My)
        Mr_ASD =  SDIComposite.C2017.Eq_A3__5b(My)

        composite_properties = SDIComposite.C2017.FlexuralProperties(

                Es,
            
                fc,
                β1,
                nothing,
            
                K1,
                K3,
                K,
                
                d,
                
                ρ,
            
                c_over_d,
                c_over_d_b,
            
                My,
                Mno_LRFD,
                Mno_ASD,
                
                nothing,
                nothing,
                nothing,
                nothing,
                    
          
            )

    elseif c_over_d >= c_over_d_b

        ϵcu = 0.003

        m = SDIComposite.C2017.EqA3__9(Es, ϵcu, fc, β1)

        c = SDIComposite.C2017.EqA3__7(d, ρ, m)

        Mro_LRFD = SDIComposite.C2017.EqA3__6a(fc, b, β1, c, d, My)

        Mro_ASD = SDIComposite.C2017.EqA3__6b(fc, b, β1, c, d, My)

        composite_properties = SDIComposite.C2017.FlexuralProperties(

        Es, 

        K1,
        K3,
        K,

        fc, 
        β1, 
        ϵcu, 
    
        d, 
    
        ρ, 

        c_over_d, 
        c_over_d_b, 

        My,
        nothing,
        nothing,
    
        m,
        c,
        Mro_LRFD,
        Mro_ASD)

    end

end




