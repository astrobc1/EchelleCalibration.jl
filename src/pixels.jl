export get_modes_gauss

function get_modes_gauss(
        spec::Vector{<:Real};
        xrange::Union{Vector{<:Real}, Nothing}=nothing,
        min_mode_spacing::Int,
        σ_bounds::Vector{<:Real},
        fit_window::Union{Vector{<:Real}, Nothing}=nothing,
        background_poly_deg::Union{Int, Nothing}=nothing,
    )
    modes0 = find_modes_centroid(spec; xrange, min_mode_spacing)
    result = fit_modes_gauss(spec, modes0; σ_bounds, background_poly_deg, fit_window)
    return result

end