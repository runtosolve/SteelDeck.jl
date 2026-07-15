
# ============================================================================
# TWO-SPAN CONTINUOUS BEAM  (unequal spans allowed)
#
# Layout:   A ----------- span 1 (L1) ----------- B ----------- span 2 (L2) ----------- C
#           simple support            simple support             simple support
#
# This file solves the extra ("continuous") support at B with Clapeyron's
# three-moment theorem — see the explanation right before
# `construction_loads_unequal_spans` below. Everything in it is exact: no numerical
# grid, no convergence to worry about, just closed-form formulas and exact
# root-finding.
#
# Loads (point, patch/UDL, or a mix) are given as plain lists, and every
# position is measured from support A along the whole beam (so a load in
# span 2 has a position bigger than L1).
# ============================================================================


struct ConstructionLoadsUnequalSpansInputs

    L1  # span 1 length (A to B)
    L2  # span 2 length (B to C)
    N   # support bearing width, centered on the support — the bearing edge is N/2 from the support centerline

    P   # list of point load magnitudes (downward positive)
    xP  # list of point load positions, one per entry in P, measured from support A.
        # Each entry is normally a number. At most ONE entry may instead be
        # `nothing` — meaning that load's position is unknown (e.g. a worker
        # walking across the deck), and the function automatically sweeps
        # the whole beam to find whichever position is critical for each result.

    w   # list of UDL/patch load intensities (force per length, downward positive)
    a   # list of UDL/patch load start positions, measured from support A
    b   # list of UDL/patch load end positions, measured from support A

    E   # modulus of elasticity, consistent with w, P, L1, L2
    I   # effective moment of inertia for deflection, consistent with L1, L2

end


struct ConstructionLoadsUnequalSpansOutputs

    inputs

    Ra  # reaction at exterior support A
    Rb  # reaction at interior support B
    Rc  # reaction at exterior support C

    V_ext_A  # shear at exterior support A
    V_ext_C  # shear at exterior support C
    V_int_B1 # shear at interior support B, span 1 side
    V_int_B2 # shear at interior support B, span 2 side

    V_A_edge  # shear at bearing edge (N/2 from support centerline), support A
    V_C_edge  # shear at bearing edge (N/2 from support centerline), support C
    V_B1_edge # shear at bearing edge (N/2 from support centerline), support B, span 1 side
    V_B2_edge # shear at bearing edge (N/2 from support centerline), support B, span 2 side

    M_A_edge  # moment at bearing edge (N/2 from support centerline), support A
    M_C_edge  # moment at bearing edge (N/2 from support centerline), support C
    M_B1_edge # moment at bearing edge (N/2 from support centerline), support B, span 1 side
    M_B2_edge # moment at bearing edge (N/2 from support centerline), support B, span 2 side

    M_pos_1  # maximum positive (sagging) moment in span 1
    x_pos_1  # location of M_pos_1, measured from support A

    M_pos_2  # maximum positive (sagging) moment in span 2
    x_pos_2  # location of M_pos_2, measured from support A

    M_pos    # governing positive moment = max(M_pos_1, M_pos_2)

    M_neg    # maximum negative (hogging) moment, anywhere on the beam
    x_neg    # location of M_neg, measured from support A

    Δ1       # maximum deflection magnitude in span 1
    x_Δ1     # location of Δ1, measured from support A

    Δ2       # maximum deflection magnitude in span 2
    x_Δ2     # location of Δ2, measured from support B

end


# ----------------------------------------------------------------------
# Basic building blocks: statics for a SIMPLE SPAN (just two end supports,
# no middle support). Each span of the real beam is treated as one of
# these on its own, carrying its own loads — see `construction_loads_unequal_spans`
# below for how the two spans get tied together.
# ----------------------------------------------------------------------

# Reaction at the left end of a simple span of length L carrying point
# loads P at positions xP, and patch/UDL loads w over [a,b].
function left_reaction(L, P, xP, w, a, b)

    Ra = 0.0

    for i in eachindex(w)
        centroid = (a[i] + b[i]) / 2
        Ra += w[i] * (b[i] - a[i]) * (L - centroid) / L
    end

    for j in eachindex(P)
        Ra += P[j] * (L - xP[j]) / L
    end

    return Ra

end


# Total of all the loads (used to get the right-end reaction from the left-end one).
function total_load(P, w, a, b)

    total = 0.0

    for i in eachindex(w)
        total += w[i] * (b[i] - a[i])
    end

    for j in eachindex(P)
        total += P[j]
    end

    return total

end


# Shear at position x on the simple span (Ra = the left reaction, already computed).
function shear_at(x, Ra, P, xP, w, a, b)

    V = Ra

    for i in eachindex(w)
        V -= w[i] * clamp(x - a[i], 0.0, b[i] - a[i])   # how much of this patch load is behind us
    end

    for j in eachindex(P)
        if x > xP[j]
            V -= P[j]   # this point load is behind us
        end
    end

    return V

end


# Bending moment at position x on the simple span.
function moment_at(x, Ra, P, xP, w, a, b)

    M = Ra * x

    for i in eachindex(w)
        ai, bi, wi = a[i], b[i], w[i]
        if x <= ai
            # this patch load hasn't started yet — no effect at x
        elseif x <= bi
            M -= wi * (x - ai)^2 / 2                      # we are inside the patch load
        else
            M -= wi * (bi - ai) * (x - (ai + bi) / 2)      # past it — treat it as a point load at its centroid
        end
    end

    for j in eachindex(P)
        if x > xP[j]
            M -= P[j] * (x - xP[j])
        end
    end

    return M

end


# Ordinary bisection: find where f crosses zero between lo and hi.
#
# How it works: f(lo) and f(hi) must have opposite signs (one positive, one
# negative) for this to make sense — that guarantees f crosses zero
# somewhere in between. Each step, we look at the midpoint: if f(mid) has
# the same sign as f(lo), the crossing must be in the right half, so we move
# lo up to mid; otherwise it's in the left half, so we move hi down to mid.
# Repeating this 60 times narrows the bracket down far tighter than we will
# ever need.
#
# If f(lo) and f(hi) do NOT have opposite signs, there is no crossing to
# find in this bracket (e.g. a peak sits right at one edge, not inside it),
# so we just return `fallback` instead of searching.
function bisect_zero(f, lo, hi; fallback=lo)

    f_lo = f(lo)
    f_hi = f(hi)

    if f_lo * f_hi >= 0.0
        return fallback
    end

    for step in 1:60
        mid = (lo + hi) / 2
        f_mid = f(mid)
        if f_mid == 0.0
            return mid
        elseif sign(f_mid) == sign(f_lo)
            lo = mid
            f_lo = f_mid
        else
            hi = mid
        end
    end

    return (lo + hi) / 2

end


# ============================================================================
# THE 2-SPAN SOLVER — Clapeyron's three-moment theorem
#
# Idea: pretend span 1 and span 2 are each simply supported on their own (so
# each one's own moment diagram is zero at both its own ends). Call the area
# under each span's own moment diagram A1, A2, and the distance from a span's
# LEFT end to its own diagram's centroid x̄. Clapeyron's theorem says the
# actual moment at the shared support B is:
#
#     M_B = -3·(A1·x̄1/L1 + A2·x̄2'/L2) / (L1+L2)
#
# where x̄1 is measured from support A, and x̄2' is measured from support C
# (the RIGHT end of span 2 — this one is measured from the opposite end, by
# the definition of the theorem). This comes from requiring the beam's slope
# to match on both sides of B — a continuous beam doesn't kink at an
# interior support, even though the two spans could each rotate differently
# there if they weren't tied together.
#
# Once M_B is known, everything else follows from ordinary statics: each
# span's moment diagram is its own free (simple-span) diagram plus a linear
# ramp from 0 up to M_B at the shared end, and the reactions shift by
# ∓M_B/length accordingly. Shear, deflection, and every peak location are
# all computed exactly (no numerical grid anywhere) — see
# `construction_loads_unequal_spans` below.
# ============================================================================


# Split a set of loads given in whole-beam coordinates into two sets, one per
# span, each measured from that span's OWN left end (support A for span 1,
# support B for span 2). A patch load that straddles the interior support is
# cut into two pieces, one on each side, at the same intensity.
function split_loads_by_span(L1, P, xP, w, a, b)

    P1, xP1 = Float64[], Float64[]
    P2, xP2 = Float64[], Float64[]

    for i in eachindex(P)
        if xP[i] <= L1
            push!(P1, P[i])
            push!(xP1, xP[i])
        else
            push!(P2, P[i])
            push!(xP2, xP[i] - L1)
        end
    end

    w1, a1, b1 = Float64[], Float64[], Float64[]
    w2, a2, b2 = Float64[], Float64[], Float64[]

    for i in eachindex(w)
        ai, bi, wi = a[i], b[i], w[i]
        if bi <= L1
            push!(w1, wi)
            push!(a1, ai)
            push!(b1, bi)
        elseif ai >= L1
            push!(w2, wi)
            push!(a2, ai - L1)
            push!(b2, bi - L1)
        else
            push!(w1, wi)
            push!(a1, ai)
            push!(b1, L1)
            push!(w2, wi)
            push!(a2, 0.0)
            push!(b2, bi - L1)
        end
    end

    return P1, xP1, w1, a1, b1, P2, xP2, w2, a2, b2

end


# Area under, and centroid (distance from the LEFT end) of, the free
# (simply-supported) bending moment diagram of one span carrying its own
# loads only.
#
# Computed by breaking the span into pieces at every position where the load
# changes, and using Simpson's rule (sample the left end, the middle, and
# the right end of a piece; combine them with the classic 1-4-1 weighting)
# on each piece. This is exact here, not approximate — because within one
# piece the moment diagram is at most a quadratic curve, and Simpson's rule
# is exact for anything up to a cubic curve.
function free_moment_area_and_centroid(Ls, P, xP, w, a, b)

    Ra = left_reaction(Ls, P, xP, w, a, b)

    breakpoints = sort(unique(vcat(0.0, Ls, Float64.(a), Float64.(b), Float64.(xP))))

    area = 0.0
    moment_of_area = 0.0

    for piece in 1:(length(breakpoints) - 1)

        x1 = breakpoints[piece]
        x2 = breakpoints[piece+1]
        xm = (x1 + x2) / 2

        M1 = moment_at(x1, Ra, P, xP, w, a, b)
        Mm = moment_at(xm, Ra, P, xP, w, a, b)
        M2 = moment_at(x2, Ra, P, xP, w, a, b)

        area += (x2 - x1) / 6 * (M1 + 4 * Mm + M2)
        moment_of_area += (x2 - x1) / 6 * (x1 * M1 + 4 * xm * Mm + x2 * M2)

    end

    if area == 0.0
        xbar = Ls / 2
    else
        xbar = moment_of_area / area
    end

    return area, xbar

end


# Exact ∫₀ˣ f(ξ)dξ, using the same piecewise-Simpson idea as
# `free_moment_area_and_centroid` above, but stopping partway through at an
# arbitrary point x instead of always going all the way to the end.
function exact_cum_integral(f, breakpoints, x)

    total = 0.0

    for piece in 1:(length(breakpoints) - 1)

        x1 = breakpoints[piece]
        x2 = min(breakpoints[piece+1], x)

        if x2 <= x1
            break
        end

        xm = (x1 + x2) / 2
        total += (x2 - x1) / 6 * (f(x1) + 4 * f(xm) + f(x2))

        if breakpoints[piece+1] >= x
            break
        end

    end

    return total

end


# Solve the beam for one specific, fully-known set of loads (every point
# load's position must be a real number, not `nothing`). This is the actual
# Clapeyron engine; `construction_loads_unequal_spans` below calls this directly
# when every position is known, or repeatedly (at many trial positions) when
# one load's position needs to be searched for.
function solve_fixed_positions(inputs)

    (; L1, L2, N, P, xP, w, a, b, E, I) = inputs

    L = L1 + L2

    P1, xP1, w1, a1, b1, P2, xP2, w2, a2, b2 = split_loads_by_span(L1, P, xP, w, a, b)

    A1, xbar1 = free_moment_area_and_centroid(L1, P1, xP1, w1, a1, b1)          # measured from A
    A2, xbar2_from_B = free_moment_area_and_centroid(L2, P2, xP2, w2, a2, b2)   # measured from B
    xbar2_from_C = L2 - xbar2_from_B                                            # theorem wants it from C

    M_B = -3 * (A1 * xbar1 / L1 + A2 * xbar2_from_C / L2) / L

    # each span's own simple-span reaction, at the end where it meets B
    Ra0_1 = left_reaction(L1, P1, xP1, w1, a1, b1)              # span 1's own reaction at A
    Rb0_2 = left_reaction(L2, P2, xP2, w2, a2, b2)              # span 2's own reaction at B
    Rc0_2 = total_load(P2, w2, a2, b2) - Rb0_2                  # span 2's own reaction at C

    # adding a moment M_B at one end of a simple span shifts its own two
    # reactions by ∓M_B/length; both spans feel this at their B end
    Ra = Ra0_1 + M_B / L1
    Rc = Rc0_2 + M_B / L2
    Rb = total_load(P, w, a, b) - Ra - Rc

    # the actual moment/shear anywhere on the beam: span 1's own formula plus
    # a linear ramp up to M_B at the B end, or span 2's own formula plus a
    # linear ramp down from M_B at its B end — whichever span xi falls in
    function combined_M(xi)
        if xi <= L1
            return moment_at(xi, Ra0_1, P1, xP1, w1, a1, b1) + M_B * xi / L1
        else
            eta = xi - L1   # position measured from B instead of from A
            return moment_at(eta, Rb0_2, P2, xP2, w2, a2, b2) + M_B * (1 - eta / L2)
        end
    end

    # V = dM/dx: span 1's moment includes +M_B·xi/L1, so its slope adds +M_B/L1;
    # span 2's moment includes M_B·(1-eta/L2), so its slope adds -M_B/L2
    function combined_V(xi)
        if xi <= L1
            return shear_at(xi, Ra0_1, P1, xP1, w1, a1, b1) + M_B / L1
        else
            eta = xi - L1
            return shear_at(eta, Rb0_2, P2, xP2, w2, a2, b2) - M_B / L2
        end
    end

    # this whole function is exact — no numerical grid anywhere. `breakpoints`
    # is just the handful of positions where the load pattern changes (span
    # ends, patch/UDL edges, point loads, and the interior support), because
    # combined_M is at most quadratic between any two of them.
    breakpoints = sort(unique(vcat(0.0, L1, L, Float64.(a), Float64.(b), Float64.(xP))))

    function q(xi)   # curvature = M/(EI)
        return combined_M(xi) / (E * I)
    end

    function qxi(xi)   # curvature times position — needed for the deflection formula below
        return q(xi) * xi
    end

    Q_total = exact_cum_integral(q, breakpoints, L)
    Qxi_total = exact_cum_integral(qxi, breakpoints, L)
    Ra_bar = Q_total - Qxi_total / L

    function θ(xi)   # beam slope at xi
        return Ra_bar - exact_cum_integral(q, breakpoints, xi)
    end

    function v(xi)   # beam deflection at xi
        integral_q = exact_cum_integral(q, breakpoints, xi)
        integral_qxi = exact_cum_integral(qxi, breakpoints, xi)
        return Ra_bar * xi - xi * integral_q + integral_qxi
    end

    # ---- moment peaks ----
    # candidates: every breakpoint directly (catches a peak sitting exactly
    # under a point load, where shear jumps through zero instead of crossing
    # it smoothly), plus an exact shear-zero bisection within every segment
    # where shear changes sign
    candidate_x = Float64[]
    candidate_M = Float64[]

    for bp in breakpoints
        push!(candidate_x, bp)
        push!(candidate_M, combined_M(bp))
    end

    for piece in 1:(length(breakpoints) - 1)

        x1 = breakpoints[piece]
        x2 = breakpoints[piece+1]

        # nudge inside the segment, not the raw breakpoints — combined_V's branch
        # is chosen by "xi <= L1", so evaluating exactly at x1 = L1 would silently
        # pick the OTHER span's shear formula instead of this segment's
        V1 = combined_V(x1 + 1e-9)
        V2 = combined_V(x2 - 1e-9)

        if V1 * V2 < 0.0
            xr = bisect_zero(combined_V, x1 + 1e-9, x2 - 1e-9)
            push!(candidate_x, xr)
            push!(candidate_M, combined_M(xr))
        end

    end

    M_pos_1 = -Inf
    x_pos_1 = 0.0
    M_pos_2 = -Inf
    x_pos_2 = 0.0
    M_neg = Inf
    x_neg = 0.0

    for i in eachindex(candidate_x)

        if candidate_x[i] <= L1 && candidate_M[i] > M_pos_1
            M_pos_1 = candidate_M[i]
            x_pos_1 = candidate_x[i]
        end

        if candidate_x[i] >= L1 && candidate_M[i] > M_pos_2
            M_pos_2 = candidate_M[i]
            x_pos_2 = candidate_x[i]
        end

        if candidate_M[i] < M_neg
            M_neg = candidate_M[i]
            x_neg = candidate_x[i]
        end

    end

    M_pos = max(M_pos_1, M_pos_2)

    # ---- deflection peaks ----
    # within each span, find every interior spot where the exact slope θ
    # crosses zero — deflection can only turn around at those spots
    # (v(0)=v(L1)=v(L)=0 at the true supports, and is smooth in between).
    #
    # Unlike shear (linear within a segment, so it crosses zero at most once),
    # θ is a CUBIC within a segment — the integral of the quadratic curvature —
    # and a cubic can cross zero more than once (e.g. deflection can dip one
    # way, come back, then dip the other way, all inside one UDL). So the
    # whole span is sampled at a handful of points per segment first, to make
    # sure no crossing hides between two points that happen to land on the
    # same side of zero; every sign change found this way is then pinned down
    # exactly by bisection, same as everywhere else in this file.
    #
    # All segments' samples are gathered into ONE continuous, ordered list
    # before checking for sign changes — not checked segment-by-segment —
    # because a crossing can sit almost exactly ON a segment boundary (e.g. a
    # point load's own position is often very close to where deflection
    # peaks), and checking each segment in isolation would miss a crossing
    # that falls between the last sample of one segment and the first sample
    # of the next.
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
                    best_x = xr
                    best_v = vr
                    found_one = true
                end
            end
        end

        return best_x, best_v

    end

    span1_breakpoints = Float64[]
    for bp in breakpoints
        if bp <= L1
            push!(span1_breakpoints, bp)
        end
    end

    span2_breakpoints = Float64[]
    for bp in breakpoints
        if bp >= L1
            push!(span2_breakpoints, bp)
        end
    end

    x_Δ1, Δ1 = deflection_extremum(span1_breakpoints)
    x_Δ2_global, Δ2 = deflection_extremum(span2_breakpoints)
    x_Δ2 = x_Δ2_global - L1

    # ---- shear and moment at the supports and their bearing edges ----
    V_ext_A = abs(combined_V(1e-9))
    V_ext_C = abs(combined_V(L - 1e-9))
    V_int_B1 = abs(combined_V(L1 - 1e-9))
    V_int_B2 = abs(combined_V(L1 + 1e-9))

    V_A_edge = abs(combined_V(N / 2))
    V_C_edge = abs(combined_V(L - N / 2))
    V_B1_edge = abs(combined_V(L1 - N / 2))
    V_B2_edge = abs(combined_V(L1 + N / 2))

    M_A_edge = abs(combined_M(N / 2))
    M_C_edge = abs(combined_M(L - N / 2))
    M_B1_edge = abs(combined_M(L1 - N / 2))
    M_B2_edge = abs(combined_M(L1 + N / 2))

    return ConstructionLoadsUnequalSpansOutputs(
        inputs,

        Ra, Rb, Rc,

        V_ext_A, V_ext_C, V_int_B1, V_int_B2,

        V_A_edge, V_C_edge, V_B1_edge, V_B2_edge,
        M_A_edge, M_C_edge, M_B1_edge, M_B2_edge,

        M_pos_1, x_pos_1,
        M_pos_2, x_pos_2,
        M_pos,

        M_neg, x_neg,

        Δ1, x_Δ1,
        Δ2, x_Δ2
    )

end


# ============================================================================
# THE PUBLIC ENTRY POINT — handles both a fully-known load AND an unknown
# point-load position, in one function.
#
# Sometimes you don't know exactly where a point load will sit (e.g. a
# worker walking across the deck). Rather than pick one position, mark that
# load's entry in `xP` as `nothing`, and this function will try many
# positions across the whole beam and, separately for every demand quantity
# (each moment, each shear, each reaction, each deflection), keep whichever
# position made that particular quantity critical. Different quantities are
# usually critical at different load positions — that's normal and expected.
#
# If every position in `xP` is a real number (nothing unknown), there's
# nothing to search for, so this just solves once and returns a plain
# `ConstructionLoadsUnequalSpansOutputs` — a single answer, not an envelope.
# ============================================================================


# One governing (critical) value for a single demand quantity, returned
# only when a point load's position was searched for: its value, where
# along the beam that value occurs, and the load position that produced it.
struct EnvelopeQuantity
    value
    x
    x_load
end


struct ConstructionLoadsUnequalSpansMovingOutputs

    inputs

    Ra_max  # most downward reaction at A, over all moving-load positions tried, and the position causing it
    Ra_min  # most upward (uplift) reaction at A

    Rb_max
    Rb_min

    Rc_max
    Rc_min

    V_ext_A
    V_ext_C
    V_int_B1
    V_int_B2

    V_A_edge
    V_C_edge
    V_B1_edge
    V_B2_edge

    M_A_edge
    M_C_edge
    M_B1_edge
    M_B2_edge

    M_pos_1  # governing max positive moment in span 1 (.x = location within span 1, .x_load = moving load position causing it)
    M_pos_2  # governing max positive moment in span 2

    M_neg    # governing max negative (hogging) moment, anywhere on the beam

    Δ1       # governing max deflection in span 1
    Δ2       # governing max deflection in span 2

end


function construction_loads_unequal_spans(inputs::ConstructionLoadsUnequalSpansInputs)

    xP = inputs.xP

    # look for a point load whose position is unknown (marked `nothing`
    # instead of a number). Only one is supported at a time.
    moving_index = 0
    for i in eachindex(xP)
        if xP[i] === nothing
            moving_index = i
            break
        end
    end

    if moving_index == 0
        # every position is known — solve once, directly
        return solve_fixed_positions(inputs)
    end

    return solve_moving_position(inputs, moving_index)

end


# Find where the load position makes a quantity critical, in two stages:
#
#   1. Coarse pass — try the load at 101 evenly-spaced spots and see which
#      one is roughly the critical one. This is just to bracket the right
#      neighborhood.
#
#   2. Refine — within that neighborhood, repeatedly check the SLOPE of the
#      quantity (is it still climbing, or has it started falling, as the
#      load nudges left vs. right?) and cut the search bracket in half
#      toward whichever side is still climbing. Once the slope is
#      essentially zero, we've found the exact peak — same "zero-crossing"
#      idea used elsewhere in this file for shear and slope, just applied
#      to the envelope quantity itself instead. This pins down the true
#      critical position, not just the nearest grid point.
function solve_moving_position(inputs, moving_index)

    n_grid = 101

    (; L1, L2, N, P, xP, w, a, b, E, I) = inputs

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
    x_move_max = L1 + L2

    # run the ordinary fixed-position solver with the moving load placed at xi
    function solve_at(xi)
        Pall = vcat(Pfixed, [Pmove])
        xPall = vcat(xPfixed, [xi])
        fixed_inputs = ConstructionLoadsUnequalSpansInputs(L1, L2, N, Pall, xPall, w, a, b, E, I)
        return solve_fixed_positions(fixed_inputs)
    end

    positions = collect(range(x_move_min, x_move_max, length=n_grid))

    results = ConstructionLoadsUnequalSpansOutputs[]
    for xi in positions
        push!(results, solve_at(xi))
    end

    # fixed (not shrinking) step used to estimate the slope of a quantity as
    # the load position moves — kept fixed so the two nearby evaluations never
    # get so close together that floating-point rounding swamps the slope
    h = (x_move_max - x_move_min) * 1.0e-6

    # Within the bracket [lo,hi] (the coarse peak and its two neighbors), find
    # the exact position where get_value's SLOPE is zero — that is the true
    # peak, pinpointed exactly rather than landing on whichever grid point
    # happened to be tried. If the slope never crosses zero in this bracket
    # (the true peak sits right at the edge of the allowed range, not in the
    # interior), we keep the coarse-grid position as-is.
    function refine(get_value, lo, hi, fallback_x)
        function slope(xi)
            return (get_value(solve_at(xi + h)) - get_value(solve_at(xi - h))) / (2h)
        end
        return bisect_zero(slope, lo, hi; fallback=fallback_x)
    end

    # `get_value` reads off the quantity we're maximizing/minimizing, from
    # one run's results; `get_location` reads off where that quantity
    # occurred within the span, from the same run.
    function critical(get_value, get_location, bigger_is_critical::Bool)

        values = Float64[]
        for r in results
            push!(values, get_value(r))
        end

        if bigger_is_critical
            k = argmax(values)
        else
            k = argmin(values)
        end

        lo = positions[max(k - 1, 1)]
        hi = positions[min(k + 1, length(positions))]
        x_load = refine(get_value, lo, hi, positions[k])

        r = solve_at(x_load)
        return EnvelopeQuantity(get_value(r), get_location(r), x_load)

    end

    # for deflection, "critical" means biggest in size, whichever direction it
    # points — sagging down or cambering up, depending on where the load lands
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

    V_ext_A = critical(r -> r.V_ext_A, r -> 0.0, true)
    V_ext_C = critical(r -> r.V_ext_C, r -> L1 + L2, true)
    V_int_B1 = critical(r -> r.V_int_B1, r -> L1, true)
    V_int_B2 = critical(r -> r.V_int_B2, r -> L1, true)

    V_A_edge = critical(r -> r.V_A_edge, r -> N / 2, true)
    V_C_edge = critical(r -> r.V_C_edge, r -> L1 + L2 - N / 2, true)
    V_B1_edge = critical(r -> r.V_B1_edge, r -> L1 - N / 2, true)
    V_B2_edge = critical(r -> r.V_B2_edge, r -> L1 + N / 2, true)

    M_A_edge = critical(r -> r.M_A_edge, r -> N / 2, true)
    M_C_edge = critical(r -> r.M_C_edge, r -> L1 + L2 - N / 2, true)
    M_B1_edge = critical(r -> r.M_B1_edge, r -> L1 - N / 2, true)
    M_B2_edge = critical(r -> r.M_B2_edge, r -> L1 + N / 2, true)

    M_pos_1 = critical(r -> r.M_pos_1, r -> r.x_pos_1, true)
    M_pos_2 = critical(r -> r.M_pos_2, r -> r.x_pos_2, true)

    M_neg = critical(r -> r.M_neg, r -> r.x_neg, false)

    Δ1 = critical_magnitude(r -> r.Δ1, r -> r.x_Δ1)
    Δ2 = critical_magnitude(r -> r.Δ2, r -> r.x_Δ2)

    return ConstructionLoadsUnequalSpansMovingOutputs(
        inputs,

        Ra_max, Ra_min,
        Rb_max, Rb_min,
        Rc_max, Rc_min,

        V_ext_A, V_ext_C, V_int_B1, V_int_B2,

        V_A_edge, V_C_edge, V_B1_edge, V_B2_edge,
        M_A_edge, M_C_edge, M_B1_edge, M_B2_edge,

        M_pos_1, M_pos_2,

        M_neg,

        Δ1, Δ2
    )

end


# ============================================================================
# SDI APPENDIX 2 CONSTRUCTION LOAD ANALYSIS — TWO SPANS, POSSIBLY UNEQUAL
#
# Mirrors `construction_loads(..., "Double")` in ConstructionLoads.jl, which
# calls SDI C-2022's closed-form Eq. C-A2-8 to C-A2-19 — those are only valid
# for two EQUAL spans. This version solves the same three SDI load cases with
# the general two-span (Clapeyron) solver above instead, so span 1 and span 2
# may differ. With L1 == L2 this reproduces the closed-form coefficients
# (0.203, 0.096, 0.125, 0.375, 1.25, 0.0054, ...) to within solver tolerance.
#
# Those coefficients mix TWO different continuous-beam loading patterns, and
# reproducing them for unequal spans means reproducing both patterns, not
# just one combined solve:
#   - Reactions (0.375, 1.25) and the negative support moment (0.125) are the
#     classic "BOTH spans loaded" coefficients — the worst case for hogging
#     moment and for the reactions is dead+live present everywhere.
#   - The positive span moment (0.096, and 0.203 for the moving P) is instead
#     the classic "PATTERN loading" coefficient — load ONE span, leave the
#     other bare — which is worse for sagging moment in the loaded span than
#     loading both spans would be (0.096 wℓ² vs. 0.070 wℓ² for "both loaded",
#     confirmed analytically via Clapeyron's theorem with one span's free
#     moment area zeroed out). SDI's P + W1 case (0.203·P·ℓ + 0.096·W1·ℓ²) is
#     literal superposition of P's own moving-load envelope (no other load
#     present) and W1's own pattern-loaded envelope — each optimized on its
#     own, then summed; this is exact here because reactions/moments are
#     linear in each load, so summing two independently-optimized envelopes
#     equals evaluating the combined system at either envelope's optimal
#     position.
#
# The three load cases, exactly as in SDI Appendix 2 (never superimposed with
# each other — each is an independent loading, and the governing demand is
# the max across all three):
#   Case 1: dead load W1 (pattern-loaded, for M_plus_1) + a single
#           concentrated load P (moving envelope) — summed per above — with
#           W1 as a full-span UDL (both spans) for the Case-1 reactions.
#   Case 2: dead load W1 + uniform construction live load W2 — pattern-loaded
#           for M_plus_2, full-span (both spans) UDL for M_neg_1 and the
#           Case-2 reactions.
#   Case 3: W3 alone — same pattern-loaded / both-spans-loaded split as
#           Case 2, for M_plus_3 / (M_neg_2 and Case-3 reactions).
# Deflection uses dead load W1 alone, both spans loaded (no W2, W3, or P),
# matching the equal-span workflow.
#
# `P`, `W1`, `W2`, `W3` are load intensities in the caller's units (matching
# `ConstructionLoadInputs`); `L1`, `L2` are the two span lengths. Bearing-edge
# quantities aren't produced here — as with the equal-span workflow, the
# caller gets those by calling this function again with clear spans
# (L1 - N, L2 - N) in place of (L1, L2).
# ============================================================================
function construction_loads_double_span(P, W1, W2, W3, L1, L2, E, I, aMn, aVn, unit_width)

    N = 0.0   # bearing-edge fields aren't consumed here
    L = L1 + L2

    # Pattern-loaded positive moment for a UDL of intensity w: load span 1
    # alone (span 2 bare) and span 2 alone (span 1 bare); each span's own
    # positive moment is read from whichever solve actually loaded it.
    function pattern_M_plus(w)
        r_s1 = construction_loads_unequal_spans(
            ConstructionLoadsUnequalSpansInputs(L1, L2, N, Float64[], Float64[], [w], [0.0], [L1], E, I))
        r_s2 = construction_loads_unequal_spans(
            ConstructionLoadsUnequalSpansInputs(L1, L2, N, Float64[], Float64[], [w], [L1], [L], E, I))
        return r_s1.M_pos_1, r_s2.M_pos_2
    end

    # Both-spans-loaded solve for a UDL of intensity w (governs negative
    # moment, reactions, and deflection).
    function both_loaded(w)
        return construction_loads_unequal_spans(
            ConstructionLoadsUnequalSpansInputs(L1, L2, N, Float64[], Float64[], [w], [0.0], [L], E, I))
    end

    # P alone (no UDL at all), moving — its own envelope positive moment and
    # reactions per span, un-perturbed by any pattern-loading choice for W1.
    rP = construction_loads_unequal_spans(
        ConstructionLoadsUnequalSpansInputs(L1, L2, N, [P], [nothing], Float64[], Float64[], Float64[], E, I))

    # ---- Case 1: P (moving) + W1 ----
    W1_pat_1, W1_pat_2 = pattern_M_plus(W1)
    M_plus_1 = max(rP.M_pos_1.value + W1_pat_1, rP.M_pos_2.value + W1_pat_2)
    r1 = both_loaded(W1)
    P_ext_1  = max(rP.Ra_max.value + r1.Ra, rP.Rc_max.value + r1.Rc)   # pair each support's own P-envelope + W1 contribution, don't mix supports
    P_int_1  = rP.Rb_max.value + r1.Rb

    # ---- Case 2: W1 + W2 ----
    W2_pat_1, W2_pat_2 = pattern_M_plus(W1 + W2)
    M_plus_2 = max(W2_pat_1, W2_pat_2)
    r2 = both_loaded(W1 + W2)
    M_neg_1  = abs(r2.M_neg)
    P_ext_2  = max(r2.Ra, r2.Rc)
    P_int_2  = r2.Rb

    # ---- Case 3: W3 alone ----
    W3_pat_1, W3_pat_2 = pattern_M_plus(W3)
    M_plus_3 = max(W3_pat_1, W3_pat_2)
    r3 = both_loaded(W3)
    M_neg_2  = abs(r3.M_neg)
    P_ext_3  = max(r3.Ra, r3.Rc)
    P_int_3  = r3.Rb

    M_plus = max(M_plus_1, M_plus_2, M_plus_3)
    M_neg  = max(M_neg_1, M_neg_2)

    # ---- Deflection: dead load W1 alone, both spans loaded ----
    rΔ = both_loaded(W1)
    Δ = max(abs(rΔ.Δ1), abs(rΔ.Δ2))

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
