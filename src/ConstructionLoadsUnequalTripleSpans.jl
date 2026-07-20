
# ============================================================================
# THREE-SPAN CONTINUOUS BEAM (unequal spans)
#
# Layout:   A --span1(L1)-- B --span2(L2)-- C --span3(L3)-- D
#           simple          simple          simple          simple
#
# Generalizes the two-span solver in ConstructionLoadsUnequalSpans.jl (see the
# long comment there for the single-equation derivation) to TWO interior
# supports, B and C, by writing Clapeyron's three-moment theorem once at EACH:
#
#   at B (spans 1,2):  2(L1+L2)·M_B +      L2 ·M_C = -6·(A1·x̄1/L1 + A2·x̄2'/L2)
#   at C (spans 2,3):       L2 ·M_B + 2(L2+L3)·M_C = -6·(A2·x̄2/L2 + A3·x̄3'/L3)
#
# where, exactly as in the two-span case, x̄ᵢ is each span's own free-moment-
# diagram centroid measured from its OWN left end, and x̄ᵢ' is the same span's
# centroid measured from its OWN right end (= Lᵢ - x̄ᵢ) — used whenever that
# span plays the "right-hand" role in an equation. Span 2 plays BOTH roles
# (right-hand span at B, left-hand span at C), so its area/centroid is
# computed once and reused both ways. This is a plain 2×2 linear system for
# (M_B, M_C); for N spans it generalizes to an (N-1)×(N-1) TRIDIAGONAL system,
# one equation per interior support, but 2 unknowns is small enough to just
# solve directly with Cramer's rule.
#
# With both M_B and M_C known, span 2 (the only span with a nonzero moment at
# BOTH ends) needs its ramp/shear generalized from "0 up to one end" (as in
# the two-span case) to a straight-line interpolation between its two end
# moments; spans 1 and 3 still have one zero end each, so their formulas look
# just like the two-span case's outer spans.
#
# `left_reaction`, `total_load`, `shear_at`, `moment_at`, `bisect_zero`,
# `free_moment_area_and_centroid`, `exact_cum_integral`, and `EnvelopeQuantity`
# are all generic (they only ever operate on ONE span's own local loads) and
# are reused as-is from ConstructionLoadsUnequalSpans.jl.
# ============================================================================


struct ConstructionLoadsUnequalTripleSpansInputs

    L1  # span 1 length (A to B)
    L2  # span 2 length (B to C)
    L3  # span 3 length (C to D)
    N   # support bearing width, centered on the support — the bearing edge is N/2 from the support centerline

    P   # list of point load magnitudes (downward positive)
    xP  # list of point load positions, one per entry in P, measured from support A.
        # Each entry is normally a number. At most ONE entry may instead be
        # `nothing` — same moving-load convention as the two-span solver.

    w   # list of UDL/patch load intensities (force per length, downward positive)
    a   # list of UDL/patch load start positions, measured from support A
    b   # list of UDL/patch load end positions, measured from support A

    E   # modulus of elasticity, consistent with w, P, L1, L2, L3
    I   # effective moment of inertia for deflection, consistent with L1, L2, L3

end


struct ConstructionLoadsUnequalTripleSpansOutputs

    inputs

    Ra  # reaction at exterior support A
    Rb  # reaction at interior support B
    Rc  # reaction at interior support C
    Rd  # reaction at exterior support D

    V_ext_A  # shear at exterior support A
    V_ext_D  # shear at exterior support D
    V_int_B1 # shear at interior support B, span 1 side
    V_int_B2 # shear at interior support B, span 2 side
    V_int_C1 # shear at interior support C, span 2 side
    V_int_C2 # shear at interior support C, span 3 side

    V_A_edge  # shear at bearing edge (N/2 from support centerline), support A
    V_D_edge  # shear at bearing edge (N/2 from support centerline), support D
    V_B1_edge # shear at bearing edge (N/2 from support centerline), support B, span 1 side
    V_B2_edge # shear at bearing edge (N/2 from support centerline), support B, span 2 side
    V_C1_edge # shear at bearing edge (N/2 from support centerline), support C, span 2 side
    V_C2_edge # shear at bearing edge (N/2 from support centerline), support C, span 3 side

    M_A_edge  # moment at bearing edge (N/2 from support centerline), support A
    M_D_edge  # moment at bearing edge (N/2 from support centerline), support D
    M_B1_edge # moment at bearing edge (N/2 from support centerline), support B, span 1 side
    M_B2_edge # moment at bearing edge (N/2 from support centerline), support B, span 2 side
    M_C1_edge # moment at bearing edge (N/2 from support centerline), support C, span 2 side
    M_C2_edge # moment at bearing edge (N/2 from support centerline), support C, span 3 side

    M_pos_1  # maximum positive (sagging) moment in span 1
    x_pos_1  # location of M_pos_1, measured from support A

    M_pos_2  # maximum positive (sagging) moment in span 2
    x_pos_2  # location of M_pos_2, measured from support A

    M_pos_3  # maximum positive (sagging) moment in span 3
    x_pos_3  # location of M_pos_3, measured from support A

    M_pos    # governing positive moment = max(M_pos_1, M_pos_2, M_pos_3)

    M_neg    # maximum negative (hogging) moment, anywhere on the beam (at B, at C, or wherever governs)
    x_neg    # location of M_neg, measured from support A

    Δ1       # maximum deflection magnitude in span 1
    x_Δ1     # location of Δ1, measured from support A

    Δ2       # maximum deflection magnitude in span 2
    x_Δ2     # location of Δ2, measured from support B

    Δ3       # maximum deflection magnitude in span 3
    x_Δ3     # location of Δ3, measured from support C

end


# Split loads given in whole-beam coordinates into three sets, one per span,
# each measured from that span's OWN left end. A patch load that straddles
# either interior support (or, in principle, both) is clipped and a piece is
# added to every span it actually overlaps.
function split_loads_by_span_triple(L1, L2, L3, P, xP, w, a, b)

    bounds = (0.0, L1, L1 + L2, L1 + L2 + L3)

    P_out  = (Float64[], Float64[], Float64[])
    xP_out = (Float64[], Float64[], Float64[])
    w_out  = (Float64[], Float64[], Float64[])
    a_out  = (Float64[], Float64[], Float64[])
    b_out  = (Float64[], Float64[], Float64[])

    for i in eachindex(P)
        x = xP[i]
        span = x <= bounds[2] ? 1 : (x <= bounds[3] ? 2 : 3)
        push!(P_out[span], P[i])
        push!(xP_out[span], x - bounds[span])
    end

    for i in eachindex(w)
        ai, bi, wi = a[i], b[i], w[i]
        for span in 1:3
            lo, hi = bounds[span], bounds[span+1]
            clipped_a = max(ai, lo)
            clipped_b = min(bi, hi)
            if clipped_b > clipped_a
                push!(w_out[span], wi)
                push!(a_out[span], clipped_a - lo)
                push!(b_out[span], clipped_b - lo)
            end
        end
    end

    return P_out[1], xP_out[1], w_out[1], a_out[1], b_out[1],
           P_out[2], xP_out[2], w_out[2], a_out[2], b_out[2],
           P_out[3], xP_out[3], w_out[3], a_out[3], b_out[3]

end


# Solve the beam for one specific, fully-known set of loads. This is the
# three-span analog of `solve_fixed_positions` above.
function solve_fixed_positions_triple(inputs)

    (; L1, L2, L3, N, P, xP, w, a, b, E, I) = inputs

    L = L1 + L2 + L3

    P1, xP1, w1, a1, b1,
    P2, xP2, w2, a2, b2,
    P3, xP3, w3, a3, b3 = split_loads_by_span_triple(L1, L2, L3, P, xP, w, a, b)

    A1, xbar1_A = free_moment_area_and_centroid(L1, P1, xP1, w1, a1, b1)   # measured from A
    A2, xbar2_B = free_moment_area_and_centroid(L2, P2, xP2, w2, a2, b2)   # measured from B
    A3, xbar3_C = free_moment_area_and_centroid(L3, P3, xP3, w3, a3, b3)   # measured from C

    xbar2_C = L2 - xbar2_B   # span 2 as the RIGHT-hand span of the equation at B
    xbar3_D = L3 - xbar3_C   # span 3 as the RIGHT-hand span of the equation at C

    # 2x2 Clapeyron system for the two interior-support moments M_B, M_C
    c11 = 2 * (L1 + L2); c12 = L2
    c21 = L2;             c22 = 2 * (L2 + L3)
    rhs1 = -6 * (A1 * xbar1_A / L1 + A2 * xbar2_C / L2)
    rhs2 = -6 * (A2 * xbar2_B / L2 + A3 * xbar3_D / L3)

    det = c11 * c22 - c12 * c21
    M_B = (rhs1 * c22 - c12 * rhs2) / det
    M_C = (c11 * rhs2 - c21 * rhs1) / det

    # each span's own simple-span reactions, under its own loads only
    Ra0_1 = left_reaction(L1, P1, xP1, w1, a1, b1)
    Rb0_1 = total_load(P1, w1, a1, b1) - Ra0_1

    Ra0_2 = left_reaction(L2, P2, xP2, w2, a2, b2)
    Rb0_2 = total_load(P2, w2, a2, b2) - Ra0_2

    Ra0_3 = left_reaction(L3, P3, xP3, w3, a3, b3)
    Rb0_3 = total_load(P3, w3, a3, b3) - Ra0_3

    # adding end moments (M_left, M_right) to a simple span shifts its own two
    # reactions by ∓(M_right - M_left)/length; span 1 has (0, M_B), span 2 has
    # (M_B, M_C), span 3 has (M_C, 0)
    Ra = Ra0_1 + M_B / L1
    Rb = (Rb0_1 - M_B / L1) + (Ra0_2 + (M_C - M_B) / L2)
    Rc = (Rb0_2 - (M_C - M_B) / L2) + (Ra0_3 - M_C / L3)
    Rd = Rb0_3 + M_C / L3

    # combined moment/shear anywhere on the beam: each span's own free formula
    # plus a straight-line ramp between its two end moments
    function combined_M(xi)
        if xi <= L1
            return moment_at(xi, Ra0_1, P1, xP1, w1, a1, b1) + M_B * xi / L1
        elseif xi <= L1 + L2
            eta = xi - L1
            return moment_at(eta, Ra0_2, P2, xP2, w2, a2, b2) + M_B * (1 - eta / L2) + M_C * (eta / L2)
        else
            eta = xi - L1 - L2
            return moment_at(eta, Ra0_3, P3, xP3, w3, a3, b3) + M_C * (1 - eta / L3)
        end
    end

    function combined_V(xi)
        if xi <= L1
            return shear_at(xi, Ra0_1, P1, xP1, w1, a1, b1) + M_B / L1
        elseif xi <= L1 + L2
            eta = xi - L1
            return shear_at(eta, Ra0_2, P2, xP2, w2, a2, b2) + (M_C - M_B) / L2
        else
            eta = xi - L1 - L2
            return shear_at(eta, Ra0_3, P3, xP3, w3, a3, b3) - M_C / L3
        end
    end

    breakpoints = sort(unique(vcat(0.0, L1, L1 + L2, L, Float64.(a), Float64.(b), Float64.(xP))))

    function q(xi)   # curvature = M/(EI)
        return combined_M(xi) / (E * I)
    end

    function qxi(xi)
        return q(xi) * xi
    end

    Q_total = exact_cum_integral(q, breakpoints, L)
    Qxi_total = exact_cum_integral(qxi, breakpoints, L)
    Ra_bar = Q_total - Qxi_total / L

    function θ(xi)
        return Ra_bar - exact_cum_integral(q, breakpoints, xi)
    end

    function v(xi)
        integral_q = exact_cum_integral(q, breakpoints, xi)
        integral_qxi = exact_cum_integral(qxi, breakpoints, xi)
        return Ra_bar * xi - xi * integral_q + integral_qxi
    end

    # ---- moment peaks ----
    candidate_x = Float64[]
    candidate_M = Float64[]

    for bp in breakpoints
        push!(candidate_x, bp)
        push!(candidate_M, combined_M(bp))
    end

    for piece in 1:(length(breakpoints) - 1)

        x1 = breakpoints[piece]
        x2 = breakpoints[piece+1]

        # nudge inside the segment — combined_V's branch is chosen by
        # "xi <= L1" / "xi <= L1+L2", so evaluating exactly at a span
        # boundary would silently pick the wrong segment's shear formula
        V1 = combined_V(x1 + 1e-9)
        V2 = combined_V(x2 - 1e-9)

        if V1 * V2 < 0.0
            xr = bisect_zero(combined_V, x1 + 1e-9, x2 - 1e-9)
            push!(candidate_x, xr)
            push!(candidate_M, combined_M(xr))
        end

    end

    M_pos_1 = -Inf; x_pos_1 = 0.0
    M_pos_2 = -Inf; x_pos_2 = 0.0
    M_pos_3 = -Inf; x_pos_3 = 0.0
    M_neg = Inf; x_neg = 0.0

    for i in eachindex(candidate_x)

        xi = candidate_x[i]
        Mi = candidate_M[i]

        if xi <= L1 && Mi > M_pos_1
            M_pos_1 = Mi; x_pos_1 = xi
        end

        if xi >= L1 && xi <= L1 + L2 && Mi > M_pos_2
            M_pos_2 = Mi; x_pos_2 = xi
        end

        if xi >= L1 + L2 && Mi > M_pos_3
            M_pos_3 = Mi; x_pos_3 = xi
        end

        if Mi < M_neg
            M_neg = Mi; x_neg = xi
        end

    end

    M_pos = max(M_pos_1, M_pos_2, M_pos_3)

    # ---- deflection peaks, per span (see the two-span solver's comment on
    # why sampling + bisection is needed: θ is cubic within a segment) ----
    function deflection_extremum(span_breakpoints)

        samples = Float64[]
        for piece in 1:(length(span_breakpoints) - 1)
            x1 = span_breakpoints[piece]
            x2 = span_breakpoints[piece+1]
            for point in range(x1 + 1e-9, x2 - 1e-9, length=21)
                push!(samples, point)
            end
        end

        theta_values = Float64[]
        for point in samples
            push!(theta_values, θ(point))
        end

        best_x = (span_breakpoints[1] + span_breakpoints[end]) / 2
        best_v = v(best_x)
        found_one = false

        for k in 1:(length(samples) - 1)
            if theta_values[k] * theta_values[k+1] < 0.0
                xr = bisect_zero(θ, samples[k], samples[k+1])
                vr = v(xr)
                if !found_one || abs(vr) > abs(best_v)
                    best_x = xr; best_v = vr; found_one = true
                end
            end
        end

        return best_x, best_v

    end

    span1_breakpoints = Float64[bp for bp in breakpoints if bp <= L1]
    span2_breakpoints = Float64[bp for bp in breakpoints if bp >= L1 && bp <= L1 + L2]
    span3_breakpoints = Float64[bp for bp in breakpoints if bp >= L1 + L2]

    x_Δ1, Δ1 = deflection_extremum(span1_breakpoints)
    x_Δ2_global, Δ2 = deflection_extremum(span2_breakpoints)
    x_Δ2 = x_Δ2_global - L1
    x_Δ3_global, Δ3 = deflection_extremum(span3_breakpoints)
    x_Δ3 = x_Δ3_global - L1 - L2

    # ---- shear and moment at the supports and their bearing edges ----
    V_ext_A  = abs(combined_V(1e-9))
    V_ext_D  = abs(combined_V(L - 1e-9))
    V_int_B1 = abs(combined_V(L1 - 1e-9))
    V_int_B2 = abs(combined_V(L1 + 1e-9))
    V_int_C1 = abs(combined_V(L1 + L2 - 1e-9))
    V_int_C2 = abs(combined_V(L1 + L2 + 1e-9))

    V_A_edge  = abs(combined_V(N / 2))
    V_D_edge  = abs(combined_V(L - N / 2))
    V_B1_edge = abs(combined_V(L1 - N / 2))
    V_B2_edge = abs(combined_V(L1 + N / 2))
    V_C1_edge = abs(combined_V(L1 + L2 - N / 2))
    V_C2_edge = abs(combined_V(L1 + L2 + N / 2))

    M_A_edge  = abs(combined_M(N / 2))
    M_D_edge  = abs(combined_M(L - N / 2))
    M_B1_edge = abs(combined_M(L1 - N / 2))
    M_B2_edge = abs(combined_M(L1 + N / 2))
    M_C1_edge = abs(combined_M(L1 + L2 - N / 2))
    M_C2_edge = abs(combined_M(L1 + L2 + N / 2))

    return ConstructionLoadsUnequalTripleSpansOutputs(
        inputs,

        Ra, Rb, Rc, Rd,

        V_ext_A, V_ext_D, V_int_B1, V_int_B2, V_int_C1, V_int_C2,

        V_A_edge, V_D_edge, V_B1_edge, V_B2_edge, V_C1_edge, V_C2_edge,
        M_A_edge, M_D_edge, M_B1_edge, M_B2_edge, M_C1_edge, M_C2_edge,

        M_pos_1, x_pos_1,
        M_pos_2, x_pos_2,
        M_pos_3, x_pos_3,
        M_pos,

        M_neg, x_neg,

        Δ1, x_Δ1,
        Δ2, x_Δ2,
        Δ3, x_Δ3
    )

end


struct ConstructionLoadsUnequalTripleSpansMovingOutputs

    inputs

    Ra_max
    Ra_min

    Rb_max
    Rb_min

    Rc_max
    Rc_min

    Rd_max
    Rd_min

    V_ext_A
    V_ext_D
    V_int_B1
    V_int_B2
    V_int_C1
    V_int_C2

    V_A_edge
    V_D_edge
    V_B1_edge
    V_B2_edge
    V_C1_edge
    V_C2_edge

    M_A_edge
    M_D_edge
    M_B1_edge
    M_B2_edge
    M_C1_edge
    M_C2_edge

    M_pos_1  # governing max positive moment in span 1 (.x = location within span 1, .x_load = moving load position causing it)
    M_pos_2
    M_pos_3

    M_neg

    Δ1
    Δ2
    Δ3

end


# Public entry point — same "solve once, or sweep a moving point load"
# behavior as `construction_loads_unequal_spans` above, generalized to three
# spans.
function construction_loads_unequal_triple_spans(inputs::ConstructionLoadsUnequalTripleSpansInputs)

    xP = inputs.xP

    moving_index = 0
    for i in eachindex(xP)
        if xP[i] === nothing
            moving_index = i
            break
        end
    end

    if moving_index == 0
        return solve_fixed_positions_triple(inputs)
    end

    return solve_moving_position_triple(inputs, moving_index)

end


function solve_moving_position_triple(inputs, moving_index)

    n_grid = 101

    (; L1, L2, L3, N, P, xP, w, a, b, E, I) = inputs

    Pmove = P[moving_index]

    Pfixed = Float64[]
    xPfixed = Float64[]
    for i in eachindex(P)
        if i != moving_index
            push!(Pfixed, P[i])
            push!(xPfixed, xP[i])
        end
    end

    x_move_min = 0.0
    x_move_max = L1 + L2 + L3

    function solve_at(xi)
        Pall = vcat(Pfixed, [Pmove])
        xPall = vcat(xPfixed, [xi])
        fixed_inputs = ConstructionLoadsUnequalTripleSpansInputs(L1, L2, L3, N, Pall, xPall, w, a, b, E, I)
        return solve_fixed_positions_triple(fixed_inputs)
    end

    positions = collect(range(x_move_min, x_move_max, length=n_grid))

    results = ConstructionLoadsUnequalTripleSpansOutputs[]
    for xi in positions
        push!(results, solve_at(xi))
    end

    h = (x_move_max - x_move_min) * 1.0e-6

    function refine(get_value, lo, hi, fallback_x)
        function slope(xi)
            return (get_value(solve_at(xi + h)) - get_value(solve_at(xi - h))) / (2h)
        end
        return bisect_zero(slope, lo, hi; fallback=fallback_x)
    end

    function critical(get_value, get_location, bigger_is_critical::Bool)

        values = Float64[]
        for r in results
            push!(values, get_value(r))
        end

        k = bigger_is_critical ? argmax(values) : argmin(values)

        lo = positions[max(k - 1, 1)]
        hi = positions[min(k + 1, length(positions))]
        x_load = refine(get_value, lo, hi, positions[k])

        r = solve_at(x_load)
        return EnvelopeQuantity(get_value(r), get_location(r), x_load)

    end

    function critical_magnitude(get_value, get_location)

        values = Float64[]
        for r in results
            push!(values, abs(get_value(r)))
        end

        k = argmax(values)

        function get_abs_value(r)
            return abs(get_value(r))
        end

        lo = positions[max(k - 1, 1)]
        hi = positions[min(k + 1, length(positions))]
        x_load = refine(get_abs_value, lo, hi, positions[k])

        r = solve_at(x_load)
        return EnvelopeQuantity(get_value(r), get_location(r), x_load)

    end

    Ra_max = critical(r -> r.Ra, r -> 0.0, true)
    Ra_min = critical(r -> r.Ra, r -> 0.0, false)

    Rb_max = critical(r -> r.Rb, r -> L1, true)
    Rb_min = critical(r -> r.Rb, r -> L1, false)

    Rc_max = critical(r -> r.Rc, r -> L1 + L2, true)
    Rc_min = critical(r -> r.Rc, r -> L1 + L2, false)

    Rd_max = critical(r -> r.Rd, r -> L1 + L2 + L3, true)
    Rd_min = critical(r -> r.Rd, r -> L1 + L2 + L3, false)

    V_ext_A  = critical(r -> r.V_ext_A, r -> 0.0, true)
    V_ext_D  = critical(r -> r.V_ext_D, r -> L1 + L2 + L3, true)
    V_int_B1 = critical(r -> r.V_int_B1, r -> L1, true)
    V_int_B2 = critical(r -> r.V_int_B2, r -> L1, true)
    V_int_C1 = critical(r -> r.V_int_C1, r -> L1 + L2, true)
    V_int_C2 = critical(r -> r.V_int_C2, r -> L1 + L2, true)

    V_A_edge  = critical(r -> r.V_A_edge,  r -> N / 2, true)
    V_D_edge  = critical(r -> r.V_D_edge,  r -> L1 + L2 + L3 - N / 2, true)
    V_B1_edge = critical(r -> r.V_B1_edge, r -> L1 - N / 2, true)
    V_B2_edge = critical(r -> r.V_B2_edge, r -> L1 + N / 2, true)
    V_C1_edge = critical(r -> r.V_C1_edge, r -> L1 + L2 - N / 2, true)
    V_C2_edge = critical(r -> r.V_C2_edge, r -> L1 + L2 + N / 2, true)

    M_A_edge  = critical(r -> r.M_A_edge,  r -> N / 2, true)
    M_D_edge  = critical(r -> r.M_D_edge,  r -> L1 + L2 + L3 - N / 2, true)
    M_B1_edge = critical(r -> r.M_B1_edge, r -> L1 - N / 2, true)
    M_B2_edge = critical(r -> r.M_B2_edge, r -> L1 + N / 2, true)
    M_C1_edge = critical(r -> r.M_C1_edge, r -> L1 + L2 - N / 2, true)
    M_C2_edge = critical(r -> r.M_C2_edge, r -> L1 + L2 + N / 2, true)

    M_pos_1 = critical(r -> r.M_pos_1, r -> r.x_pos_1, true)
    M_pos_2 = critical(r -> r.M_pos_2, r -> r.x_pos_2, true)
    M_pos_3 = critical(r -> r.M_pos_3, r -> r.x_pos_3, true)

    M_neg = critical(r -> r.M_neg, r -> r.x_neg, false)

    Δ1 = critical_magnitude(r -> r.Δ1, r -> r.x_Δ1)
    Δ2 = critical_magnitude(r -> r.Δ2, r -> r.x_Δ2)
    Δ3 = critical_magnitude(r -> r.Δ3, r -> r.x_Δ3)

    return ConstructionLoadsUnequalTripleSpansMovingOutputs(
        inputs,

        Ra_max, Ra_min,
        Rb_max, Rb_min,
        Rc_max, Rc_min,
        Rd_max, Rd_min,

        V_ext_A, V_ext_D, V_int_B1, V_int_B2, V_int_C1, V_int_C2,

        V_A_edge, V_D_edge, V_B1_edge, V_B2_edge, V_C1_edge, V_C2_edge,
        M_A_edge, M_D_edge, M_B1_edge, M_B2_edge, M_C1_edge, M_C2_edge,

        M_pos_1, M_pos_2, M_pos_3,

        M_neg,

        Δ1, Δ2, Δ3
    )

end


# ============================================================================
# SDI APPENDIX 2 CONSTRUCTION LOAD ANALYSIS — THREE SPANS, POSSIBLY UNEQUAL
#
# Mirrors `construction_loads_double_span` in ConstructionLoadsUnequalSpans.jl
# (see that function's long comment for the "both/all loaded" vs "pattern
# loaded" distinction), extended to three spans and TWO interior supports.:
#
#   - REACTIONS at every support are worst with ALL THREE spans loaded — more
#     total load on the beam can only add to (never relieve) a support's
#     total reaction. Verified analytically: for equal spans ℓ, "all loaded"
#     gives M_B = M_C = -0.100 wℓ² and reactions 0.4 wℓ (exterior) / 1.1 wℓ
#     (interior) — matching SDI's Eq. C-A2-25-to-30 coefficients (0.4, 1.1)
#     exactly.
#   - NEGATIVE MOMENT at an interior support is instead worst with only the
#     TWO SPANS ADJACENT to that support loaded, and the FAR span left bare
#     (the classic ACI alternate-span pattern for negative moment) — loading
#     the far span actually relieves this particular support's own hogging
#     moment via the coupled 2×2 system. Verified analytically: "spans 1&2
#     loaded, span 3 bare" gives M_B = -0.1167 wℓ² — matching SDI's 0.117
#     coefficient (Eq. C-A2-20-to-24), NOT the "all loaded" value of 0.100.
#   - POSITIVE MOMENT in an end span is worst with BOTH end spans loaded and
#     the middle span bare (alternate-span pattern); positive moment in the
#     middle span is worst with the middle span loaded ALONE. Same
#     checkerboard principle SDI's own equal-span coefficients are built on.
#   - Deflection uses dead load W1 alone, ALL THREE spans loaded (same
#     principle as reactions — more load can only increase net downward
#     deflection under self-weight alone).
#
# The three SDI load cases (never superimposed with each other — the
# governing demand is the max across all three), exactly as in the two-span
# version:
#   Case 1: dead load W1 (patterned per above) + a single concentrated load P
#           (moving envelope, no other load present).
#   Case 2: dead load W1 + uniform construction live load W2.
#   Case 3: W3 alone.
# ============================================================================
function construction_loads_triple_span(P, W1, W2, W3, L1, L2, L3, E, I, aMn, aVn, unit_width)

    N = 0.0   # bearing-edge fields aren't consumed here
    L = L1 + L2 + L3
    bB = L1        # position of support B, whole-beam coordinates
    bC = L1 + L2   # position of support C, whole-beam coordinates

    # Positive-moment pattern for a UDL of intensity w: load ONE span alone,
    # the other two bare — verified analytically against SDI's own published
    # coefficients (NOT the true ACI checkerboard/alternate-span pattern,
    # which gives an even higher, more conservative number: for equal spans,
    # "span 1 alone" gives 0.0939 wℓ² — matching SDI's 0.094 almost exactly —
    # versus 0.10125 wℓ² for "spans 1&3 loaded, 2 bare". SDI's own tables use
    # the simpler single-span-loaded assumption, same convention as the
    # two-span case's 0.096 coefficient, so this generalizes that convention
    # rather than introducing a stricter (and inconsistent) one.)
    function pattern_M_plus(w)
        r1 = construction_loads_unequal_triple_spans(
            ConstructionLoadsUnequalTripleSpansInputs(L1, L2, L3, N, Float64[], Float64[], [w], [0.0], [bB], E, I))
        r2 = construction_loads_unequal_triple_spans(
            ConstructionLoadsUnequalTripleSpansInputs(L1, L2, L3, N, Float64[], Float64[], [w], [bB], [bC], E, I))
        r3 = construction_loads_unequal_triple_spans(
            ConstructionLoadsUnequalTripleSpansInputs(L1, L2, L3, N, Float64[], Float64[], [w], [bC], [L], E, I))
        return r1.M_pos_1, r2.M_pos_2, r3.M_pos_3
    end

    # Negative-moment pattern at support B: spans 1&2 loaded, span 3 bare.
    function pattern_M_neg_B(w)
        r = construction_loads_unequal_triple_spans(
            ConstructionLoadsUnequalTripleSpansInputs(L1, L2, L3, N, Float64[], Float64[],
                [w], [0.0], [bC], E, I))
        return r.M_neg   # global min already lands at B under this pattern
    end

    # Negative-moment pattern at support C: spans 2&3 loaded, span 1 bare.
    function pattern_M_neg_C(w)
        r = construction_loads_unequal_triple_spans(
            ConstructionLoadsUnequalTripleSpansInputs(L1, L2, L3, N, Float64[], Float64[],
                [w], [bB], [L], E, I))
        return r.M_neg   # global min already lands at C under this pattern
    end

    # All three spans loaded (governs reactions and deflection).
    function all_loaded(w)
        return construction_loads_unequal_triple_spans(
            ConstructionLoadsUnequalTripleSpansInputs(L1, L2, L3, N, Float64[], Float64[], [w], [0.0], [L], E, I))
    end

    # P alone (no UDL at all), moving — its own envelope, per span/support.
    rP = construction_loads_unequal_triple_spans(
        ConstructionLoadsUnequalTripleSpansInputs(L1, L2, L3, N, [P], [nothing], Float64[], Float64[], Float64[], E, I))

    # ---- Case 1: P (moving) + W1 ----
    w1_1, w1_2, w1_3 = pattern_M_plus(W1)
    M_plus_1 = max(
        rP.M_pos_1.value + w1_1,
        rP.M_pos_2.value + w1_2,
        rP.M_pos_3.value + w1_3,
    )
    r1 = all_loaded(W1)
    P_ext_1 = max(rP.Ra_max.value + r1.Ra, rP.Rd_max.value + r1.Rd)   # pair each support's own P-envelope + W1 contribution
    P_int_1 = max(rP.Rb_max.value + r1.Rb, rP.Rc_max.value + r1.Rc)

    # ---- Case 2: W1 + W2 ----
    w2_1, w2_2, w2_3 = pattern_M_plus(W1 + W2)
    M_plus_2 = max(w2_1, w2_2, w2_3)
    M_neg_1  = max(abs(pattern_M_neg_B(W1 + W2)), abs(pattern_M_neg_C(W1 + W2)))
    r2 = all_loaded(W1 + W2)
    P_ext_2  = max(r2.Ra, r2.Rd)
    P_int_2  = max(r2.Rb, r2.Rc)

    # ---- Case 3: W3 alone ----
    w3_1, w3_2, w3_3 = pattern_M_plus(W3)
    M_plus_3 = max(w3_1, w3_2, w3_3)
    M_neg_2  = max(abs(pattern_M_neg_B(W3)), abs(pattern_M_neg_C(W3)))
    r3 = all_loaded(W3)
    P_ext_3  = max(r3.Ra, r3.Rd)
    P_int_3  = max(r3.Rb, r3.Rc)

    M_plus = max(M_plus_1, M_plus_2, M_plus_3)
    M_neg  = max(M_neg_1, M_neg_2)

    # ---- Deflection: dead load W1 alone, all three spans loaded ----
    rΔ = all_loaded(W1)
    Δ = max(abs(rΔ.Δ1), abs(rΔ.Δ2), abs(rΔ.Δ3))

    Mbar = max(M_neg, M_plus)
    Vbar = max(P_ext_1, P_ext_2, P_ext_3, P_int_1, P_int_2, P_int_3)
    interaction = AISIS100.v16S3.h21(Mbar / unit_width, Vbar / unit_width, aMn, aVn)

    return (
        M_plus_1 = M_plus_1, M_plus_2 = M_plus_2, M_plus_3 = M_plus_3, M_plus = M_plus,
        M_neg_1  = M_neg_1,  M_neg_2  = M_neg_2,  M_neg  = M_neg,
        P_ext_1  = P_ext_1,  P_ext_2  = P_ext_2,  P_ext_3  = P_ext_3,
        P_int_1  = P_int_1,  P_int_2  = P_int_2,  P_int_3  = P_int_3,
        Δ = Δ,
        interaction = interaction,
    )

end
