export fit_modes_gauss

function fit_mode_gauss(
        spec::Vector{<:Real}, mode0::Real, fit_window::Vector{<:Real};
        μ_bounds::Vector{<:Real},
        σ_bounds::Vector{<:Real},
        background_poly_deg::Union{Int, Nothing}=nothing
    )

    # Window
    nx = length(spec)
    xx = max(round(mode0 + fit_window[1]), 1):min(round(mode0 + fit_window[2]), nx)
    yy = spec[Int.(xx)]

    # Remove approx baseline and normalize
    if !isnothing(background_poly_deg)
        background = nanminimum(yy)
        yy .-= background
    else
        background = nothing
    end
    scale = nanmaximum(yy)
    yy ./= scale

    # Initial parameters
    # Amp, mean, sigma, c0, c1, c2, ...
    p0 = [1.0, mode0, nanmean(σ_bounds)]
    lb = [0.5, mode0 + μ_bounds[1], σ_bounds[1]]
    ub = [1.5, mode0 + μ_bounds[2], σ_bounds[2]]

    if !isnothing(background_poly_deg)
        push!(p0, 0.001)
        push!(lb, -0.5)
        push!(ub, 0.5)
        for i=1:background_poly_deg
            push!(p0, 0.001)
            push!(lb, -0.5)
            push!(ub, 0.5)
        end
    end

    # Model
    if !isnothing(background_poly_deg)
        model = (x, p) -> begin
            return p[1] .* exp.(-0.5 .* ((x .- p[2]) ./ p[3]).^2) .+ Polynomial(p[4:end]).(x .- nanmean(x))
        end
    else
        model = (x, p) -> begin
            return p[1] .* exp.(-0.5 .* ((x .- p[2]) ./ p[3]).^2)
        end
    end

    # Good data
    good = findall(isfinite.(yy))
    if length(good) < length(p0)
        println("Could not fit mode at x=$(mode0), insufficient dof")
        return nothing
    end

    # Fit
    result = LsqFit.curve_fit(model, xx[good], yy[good], p0, lower=lb, upper=ub, autodiff=:forwarddiff)
    
    # RMS
    rms = sqrt(nansum(result.resid.^2) / length(result.resid))

    # Errors
    pbest = result.param
    pbest_err = nothing
    try
        pbest_err = LsqFit.standard_errors(result)
    catch
        println("Could not fit mode at x=$(mode0), failed errors (bad jacobian)")
        return nothing
    end

    # Pack results
    out = (;
        amp=pbest[1] * scale, amp_err=pbest_err[1] * scale,
        μ=pbest[2], μ_err=pbest_err[2],
        σ=pbest[3], σ_err=pbest_err[3],
        background_coeffs=pbest[4:end], background_coeffs_err=pbest_err[4:end],
        rms=rms * scale,
        scale,
        background
    )

    # Return
    return out

end

function fit_modes_gauss(
        spec::Vector{<:Real},
        modes::Vector{<:Real};
        background_poly_deg::Union{Int, Nothing}=nothing,
        fit_window::Union{Vector{<:Real}, Nothing}=nothing,
        σ_bounds::Vector{<:Real}=[0.5, 10.0]
    )

    # Fit results
    n_modes = length(modes)
    amplitudes = fill(NaN, n_modes)
    amplitudes_err = fill(NaN, n_modes)
    pixels = fill(NaN, n_modes)
    pixels_err = fill(NaN, n_modes)
    σs = fill(NaN, n_modes)
    σs_err = fill(NaN, n_modes)
    backgrounds = fill(NaN, n_modes)
    background_polys = Vector{Polynomial}(undef, n_modes)
    background_polys_err = Vector{Vector{Float64}}(undef, n_modes)
    rms = fill(NaN, n_modes)

    # Copy the spectrum
    spec = copy(spec)

    ds = diff(modes)
    push!(ds, ds[end])

    # Fit peaks
    for i=1:n_modes

        # Fit window
        if !isnothing(fit_window)
            _fit_window = fit_window
        else
            _fit_window = [-ds[i] / 2.25, ds[i] / 2.25]
        end


        # Fit
        result = fit_mode_gauss(spec, modes[i], _fit_window; σ_bounds, μ_bounds=[-1, 1], background_poly_deg)

        # Store results
        if !isnothing(result)
            amplitudes[i] = result.amp
            amplitudes_err[i] = result.amp_err
            pixels[i] = result.μ
            pixels_err[i] = result.μ_err
            σs[i] = result.σ
            σs_err[i] = result.σ_err
            backgrounds[i] = result.background
            background_polys[i] = Polynomial(result.background_coeffs)
            background_polys_err[i] = result.background_coeffs_err
            rms[i] = result.rms
        end

    end

    # Return
    return (;pixels, pixels_err, amplitudes, amplitudes_err, σs, σs_err, background_polys, background_polys_err, rms, backgrounds)

end