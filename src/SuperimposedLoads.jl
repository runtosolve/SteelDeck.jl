
struct SuperimposedLoadInputs 

    E
    Id 
     
    Mn_ASD
    Vn_ASD 

    panel_thickness
    total_slab_depth
    span_length 
    
    W1_steel 
    W1_concrete 
    
    Δ_over_L_limit
    
end


struct SuperimposedLoadOutputs

    inputs

    w_u_M 
    w_u_V 
    w_u_Δ

    w_u 

    w_u_limit_state 

    deck_slab_weight 

    superimposed_load 


end




function calculate_superimposed_load(inputs)

    (;
    E,
    Id, 
     
    Mn_ASD,
    Vn_ASD, 

    panel_thickness,
    total_slab_depth,
    span_length, 
    
    W1_steel, 
    W1_concrete, 
    
    Δ_over_L_limit) = inputs  


    w_u_M = Mn_ASD * 8 / (span_length * 12).^2

    w_u_M = w_u_M .* 12 #convert to lbs/ft

    w_u_V = Vn_ASD * 2 / span_length

    w_u_Δ = (Δ_over_L_limit * span_length * 12)  * ((384 / 5) * (E * Id)) / (span_length * 12)^4  * 12  #convert to lbs/ft 

    w_u = minimum([w_u_M, w_u_V, w_u_Δ])

    w_u_index = argmin([w_u_M, w_u_V, w_u_Δ])

    limit_states = ["flexure", "shear", "deflection"]

    w_u_limit_state = limit_states[w_u_index]

    deck_slab_weight = W1_concrete .+ W1_steel

    superimposed_load = minimum([w_u - deck_slab_weight, 400.0]) 

    outputs = SuperimposedLoadOutputs(

        inputs,

        w_u_M, 
        w_u_V, 
        w_u_Δ,

        w_u, 

        w_u_limit_state, 

        deck_slab_weight, 

        superimposed_load
        ) 

    return outputs 


end