
export get_lfc_modes


function get_lfc_modes(
        λ::Vector{<:Real},
        spec::Vector{<:Real},
        ν0::Real, Δν::Real;
        σ_bounds::Vector{<:Real},
        xrange::Union{Vector{Int}, Nothing}=nothing,
        background_poly_deg::Union{Int, Nothing}=nothing,
        fit_window::Union{Vector{<:Real}, Nothing}=nothing
    )

    # Resolve xrange
    spec = copy(spec)
    spec[.~(xrange[1] .<= eachindex(spec) .<= xrange[2])] .= NaN
    xi, xf = xrange
    xarr = 1:length(spec)

    # Wave range
    λi, λf = λ[xi], λ[xf]

    # Gen lfc modes
    modes_λ, _, integers = gen_lfc_modes(ν0, Δν; λi, λf)

    # Ignore first and last mode just in case
    modes_λ = modes_λ[2:end-1]
    integers = integers[2:end-1]

    # Initial modes in pixel space
    modes_pix = [Float64.(xarr[nanargmin(abs.(λ .- λi))]) for λi in modes_λ]
    ds = diff(modes_pix)
    push!(ds, ds[end])
    
    # Refine from centroids
    for i in eachindex(modes_pix)
        modes_pix[i] = refine_centroid1d(spec, modes_pix[i], [-ds[i] /  2, ds[i] / 2]; n_iterations=3)
    end

    # Fit
    results = fit_modes_gauss(spec, modes_pix; background_poly_deg, σ_bounds, fit_window)

    # Bundle results
    out = merge(results, (;λs=modes_λ, integers))

    # Return
    return out

end



function get_lfc_modes_from_known_line(
        spec::Vector{<:Real};
        x0::Real, λ0::Real,
        ν0::Real, Δν::Real,
        σ_bounds::Vector{<:Real},
        xrange::Union{Vector{Int}, Nothing}=nothing,
        background_poly_deg::Union{Int, Nothing}=nothing,
        min_mode_spacing::Real,
    )

    # Resolve xrange
    spec = copy(spec)
    spec[.~(xrange[1] .<= eachindex(spec) .<= xrange[2])] .= NaN

    # Find peaks in pixel space
    result = get_modes_gauss(spec; xrange, min_mode_spacing, σ_bounds, background_poly_deg)

    # Mode closest to this wavelength
    λs, _, integers = gen_lfc_modes(ν0, Δν; λi = λ0 - 200, λf = λ0 + 200)
    dλs = λs .- λ0
    k0 = nanargmin(abs.(dλs))

    # Label modes
    n_modes = length(result.pixels)
    dxs = x0 .- result.pixels
    j0 = nanargmin(abs.(dxs))
    ints_out = integers[(k0 - j0 + 1):(k0 + (n_modes-j0))]
    λs_out = gen_lfc_modes(ν0, Δν, ints_out).λs

    # Out
    out = merge(result, (;integers=ints_out, λs=λs_out))

    # Return
    return out

end



function gen_lfc_modes(ν0::Real, Δν::Real; λi::Real, λf::Real)
    νf = (SPEED_OF_LIGHT_MPS / λi) * 1E9
    νi = (SPEED_OF_LIGHT_MPS / λf) * 1E9
    n1 = Int(ceil((νi - ν0) / Δν))
    n2 = Int(ceil((νf - ν0) / Δν))
    integers = collect(min(n1, n2):max(n1, n2))
    λs, νs = gen_lfc_modes(ν0, Δν, integers)
    reverse!(λs)
    reverse!(νs)
    reverse!(integers)
    return (;λs, νs, integers)
end


function gen_lfc_modes(ν0::Real, Δν::Real, integers::AbstractVector{Int})
    νs = @. integers * Δν + ν0
    λs = SPEED_OF_LIGHT_MPS ./ νs .* 1E9
    return (;λs, νs)
end


function refine_centroid1d(y::Vector{<:Real}, x0::Real, window::Vector{<:Real}; n_iterations::Int=3)
    xc = x0
    nx = length(y)
    xarr = 1:nx
    for _=1:n_iterations
        k1 = max(1, Int(round(xc + window[1])))
        k2 = min(nx, Int(round(xc + window[2])))
        xx, yy = xarr[k1:k2], view(y, k1:k2)
        xc_new = nansum(xx .* yy) / nansum(yy)
        if xc_new == xc
            break
        else
            xc = xc_new
        end
    end
    return xc
end