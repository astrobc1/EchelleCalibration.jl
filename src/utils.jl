


function estimate_background(spec::Vector{<:Real}, modes::Vector{<:Real})
    nx = length(spec)
    n_modes = length(modes)
    valleys = fill(NaN, n_modes + 1)
    valleys[1] = modes[1] - (modes[2] - modes[1]) / 2
    valleys[end] = modes[end] + (modes[end] - modes[end-1]) / 2
    for i=2:(length(valleys)-1)
        valleys[i] = (modes[i-1] + modes[i]) / 2
    end
    spec_smooth = quantile_filter1d(spec, width=3)
    spec_interp_valleys = lin_interp(1:nx, spec_smooth, valleys)
    spl = CubicSpline(spec_interp_valleys, valleys)
    background = spl.(1:nx)
    return background
end

function estimate_background(spec::Vector{<:Real}, min_mode_spacing::Real)
    nx = length(spec)
    spec_smooth = quantile_filter1d(spec, width=3)
    baseline = quantile_filter1d(spec_smooth, width=nextodd(2 * min_mode_spacing), p=0)
    good = findall(isfinite.(baseline))
    background = savitzky_golay(baseline, nextodd(3 * min_mode_spacing + 1), 2).y
    return background
end

function get_mode_spacing_poly(modes::Vector{<:Real})
    d = diff(modes)
    xs = fill(NaN, 10)
    dys = fill(NaN, 10)
    bins = range(modes[1], modes[end], length=11)
    mv = modes[1:end-1]
    for i=1:10
        _xi, _xf = Int(round(bins[i])), Int(round(bins[i+1]))
        xs[i] = (_xi + _xf) / 2
        use = findall(mv .> _xi .&& mv .<= _xf)
        if length(use) > 0
            dys[i] = nanquantile(d[Int.(round.(use))], 0.25)
        end
    end
    good = findall(isfinite.(dys))
    mode_spacing_poly = @views Polynomials.fit(xs[good], dys[good], 2)
    return mode_spacing_poly
end