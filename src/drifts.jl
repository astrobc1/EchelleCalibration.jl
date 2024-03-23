
function compute_drifts_like_modes(
        modes_pixels0::Vector{<:Real}, modes_pixels_err0::Vector{<:Real},
        modes_pixels::Vector{<:Real}, modes_pixels_err::Vector{<:Real}; thresh::Real=3.0
    )
    pixels = Float64[]
    drifts = Float64[]
    for i in eachindex(modes_pixels)
        ds = modes_pixels[i] .- modes_pixels0
        k = nanargmin(abs.(ds))
        if abs(ds[k]) > thresh
            continue
        end
        push!(pixels, (modes_pixels[i] + modes_pixels0[k]) / 2)
        push!(drifts, ds[k])
        push!(drifts_err, sqrt(modes_pixels_err[i]^2 + modes_pixels_err0[k]^2))
    end
    return pixels, drifts, drifts_err
end