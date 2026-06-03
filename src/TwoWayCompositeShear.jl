struct TwoWayShearInputs

    c1
    c2

    hc

    fc

    design_method

end



struct TwoWayShearOutputs

    inputs

    βc
    bo
    Vpr
    aVpr

end



function calculate_two_way_composite_shear_strength(inputs)

    (; c1, c2, hc, fc, design_method) = inputs

    d = hc

    βc = maximum([c1, c2]) / minimum([c1, c2])

    bo = 2 * (c1 + d) + 2 * (c2 + d)

    Vpr, aVpr = SDIComposite.C2022.EqF_5_1a(βc, bo, hc, fc, design_method)

    outputs = SteelDeck.TwoWayShearOutputs(
                inputs,
                βc,
                bo,
                Vpr,
                aVpr
    )

    return outputs

end
