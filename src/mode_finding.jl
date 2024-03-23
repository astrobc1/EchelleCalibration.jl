export find_modes_centroid

function find_modes_centroid(
        spec::Vector{<:Real};
        xrange::Union{Vector{Int}, Nothing}=nothing, min_mode_spacing::Union{Real, Nothing}=nothing,
        height::Union{Real, Nothing}=nothing
    )

    # Resolve xrange
    spec = copy(spec)
    spec[.~(xrange[1] .<= eachindex(spec) .<= xrange[2])] .= NaN

    # Smooth
    spec_smooth = quantile_filter(spec, window=3)

    # Find peak indices
    peaks, _ = findpeaks1d(spec_smooth; height, distance=min_mode_spacing)

    # Ignore first and last peak
    peaks = peaks[2:end-1]

    # Iteratively refine based on centroid
    modes = refine_modes_centroid(spec_smooth, peaks, min_mode_spacing, n_iterations=3)

    # Return
    return modes

end


function refine_modes_centroid(spec, modes, window::Real; n_iterations=3)
    nx = length(spec)
    modes_out = float.(modes)
    xarr = 1:nx
    w2 = window / 2
    for _=1:n_iterations
        for j in eachindex(modes)
            use = findall((xarr .>= modes[j] - w2) .&& (xarr .<= modes[j] + w2))
            xx, yy = @views xarr[use], spec[use]
            modes_out[j] = nansum(xx .* yy) / nansum(yy)
        end
    end
    return modes_out
end

function refine_modes_centroid(spec, modes, windows::AbstractVector{<:Real}; n_iterations=3)
    nx = length(spec)
    modes_out = float.(modes)
    xarr = 1:nx
    for _=1:n_iterations
        for j in eachindex(modes)
            w2 = windows[j] / 2
            use = findall((xarr .>= modes[j] - w2) .&& (xarr .<= modes[j] + w2))
            xx, yy = @views xarr[use], spec[use]
            modes_out[j] = nansum(xx .* yy) / nansum(yy)
        end
    end
    return modes_out
end