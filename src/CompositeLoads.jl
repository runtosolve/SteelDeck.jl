
struct CompositeLoadInputs

    L   # span length
    N   # support bearing width
    w   # vector of UDL intensities (force per length)
    a   # vector of UDL start positions from left support
    b   # vector of UDL end positions from left support
    E   # modulus of elasticity, consistent with w, L (e.g. force/length^2)
    I   # effective moment of inertia for deflection, consistent with L (e.g. length^4)

end


struct CompositeLoadOutputs

    inputs

    Ra      # left support reaction
    Rb      # right support reaction
    V_max   # maximum shear (at supports)
    M_max   # maximum bending moment
    x_max_M # location of maximum bending moment
    V_edge  # shear at support edge (governing)
    M_edge  # moment at support edge (governing)
    Δ       # maximum deflection (magnitude), from the loads in `inputs.w`
    x_max_Δ # location of maximum deflection

end


function composite_loads(inputs::CompositeLoadInputs)

    (; L, N, w, a, b, E, I) = inputs

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

    V_edge = max(abs(V(N)), abs(V(L - N)))
    M_edge = max(abs(M(N)), abs(M(L - N)))

    # Maximum deflection via the conjugate-beam method: q(x) = M(x)/(E·I) is applied as a fictitious distributed load on a conjugate simply supported beam of the same span. 
    # The conjugate beam's shear equals the real beam's slope θ(x), and its moment equals the real beam's deflection v(x).
    #
    # For downward-only loads, M(x) ≥ 0 everywhere on a simple span, so θ(x) = Ra_bar - ∫₀ˣ q(ξ)dξ is monotonically non-increasing with exactly one zero crossing — the location of the single interior maximum of v(x)
    # (v is concave with v(0) = v(L) = 0).

    q(x) = M(x) / (E * I)
    qxi(x) = q(x) * x

    # exact ∫₀ˣ f(ξ)dξ via piecewise Simpson's rule; exact for any polynomial up to degree 3, which covers q (quadratic) and qxi (cubic) here since M(x) is at most quadratic within any interval bounded by x_crit
    function cum_integral(f, x)
        total = 0.0
        for j in 1:(length(x_crit) - 1)
            x1 = x_crit[j]
            x2 = min(x_crit[j+1], x)
            x2 <= x1 && break
            xm = (x1 + x2) / 2
            total += (x2 - x1) / 6 * (f(x1) + 4 * f(xm) + f(x2))
            x_crit[j+1] >= x && break
        end
        return total
    end

    Q_total   = cum_integral(q, L)
    Qxi_total = cum_integral(qxi, L)
    Ra_bar    = Q_total - Qxi_total / L

    θ(x) = Ra_bar - cum_integral(q, x)
    v(x) = Ra_bar * x - x * cum_integral(q, x) + cum_integral(qxi, x)

    x_zero_θ = Float64[]
    for j in 1:(length(x_crit) - 1)
        x1, x2 = x_crit[j], x_crit[j+1]
        θ1, θ2 = θ(x1), θ(x2)
        if θ1 == 0.0
            push!(x_zero_θ, x1)
        elseif θ2 == 0.0
            push!(x_zero_θ, x2)
        elseif θ1 * θ2 < 0.0
            lo, hi = x1, x2
            θ_lo = θ1
            for _ in 1:60
                mid = (lo + hi) / 2
                θ_mid = θ(mid)
                if θ_mid == 0.0
                    lo = hi = mid
                    break
                elseif sign(θ_mid) == sign(θ_lo)
                    lo, θ_lo = mid, θ_mid
                else
                    hi = mid
                end
            end
            push!(x_zero_θ, (lo + hi) / 2)
        end
    end

    Δ_vals  = v.(x_zero_θ)
    Δ       = isempty(Δ_vals) ? v(L / 2) : maximum(Δ_vals)
    x_max_Δ = isempty(Δ_vals) ? L / 2 : x_zero_θ[argmax(Δ_vals)]

    return CompositeLoadOutputs(inputs, Ra, Rb, V_max, M_max, x_max_M, V_edge, M_edge, Δ, x_max_Δ)

end