
struct CompositeLoadInputs

    L   # span length
    w   # vector of UDL intensities (force per length)
    a   # vector of UDL start positions from left support
    b   # vector of UDL end positions from left support

end


struct CompositeLoadOutputs

    inputs

    Ra      # left support reaction
    Rb      # right support reaction
    V_max   # maximum shear (at supports)
    M_max   # maximum bending moment
    x_max_M # location of maximum bending moment

end


function composite_loads(inputs::CompositeLoadInputs)

    (; L, w, a, b) = inputs

    n = length(w)

    Ra = sum(w[i] * (b[i] - a[i]) * (L - (a[i] + b[i]) / 2) / L for i in 1:n)
    Rb = sum(w[i] * (b[i] - a[i]) for i in 1:n) - Ra

    function V(x)
        return Ra - sum(w[i] * clamp(x - a[i], 0.0, b[i] - a[i]) for i in 1:n)
    end

    function M(x)
        function m(wi, ai, bi)
            if     x <= ai;  return 0.0
            elseif x <= bi;  return -wi * (x - ai)^2 / 2
            else;            return -wi * (bi - ai) * (x - (ai + bi) / 2)
            end
        end
        return Ra * x + sum(m(w[i], a[i], b[i]) for i in 1:n)
    end

    x_crit = sort(unique(vcat(0.0, collect(a), collect(b), L)))

    x_zero_V = Float64[]
    for j in 1:(length(x_crit) - 1)
        x1, x2 = x_crit[j], x_crit[j+1]
        V1, V2 = V(x1 + 1e-12), V(x2 - 1e-12)
        if V1 * V2 <= 0
            x_z = x1 + V1 / (V1 - V2) * (x2 - x1)
            push!(x_zero_V, clamp(x_z, x1, x2))
        end
    end

    M_vals  = M.(x_zero_V)
    M_max   = isempty(M_vals) ? 0.0 : maximum(M_vals)
    x_max_M = isempty(M_vals) ? L / 2 : x_zero_V[argmax(M_vals)]

    V_max = max(abs(Ra), abs(Rb))

    return CompositeLoadOutputs(inputs, Ra, Rb, V_max, M_max, x_max_M)

end