using PyCall

function build_λsol2d_poly_flat(coeffs, pixels, orders, min_pixel, max_pixel, min_order, max_order)
    numpy = pyimport("numpy")
    pixels_norm = normalize_val.(pixels, min_pixel, max_pixel)
    orders_norm = normalize_val.(orders, min_order, max_order)
    return numpy.polynomial.polynomial.polyval2d(pixels_norm, orders_norm, coeffs) ./ orders
end

function build_λsol2d_poly_full(coeffs, pixels, orders, min_pixel, max_pixel, min_order, max_order)
    out = Vector{Vector{Float64}}(undef, length(orders))
    for (i, order) in enumerate(orders)
        out[i] = build_λsol2d_poly_flat(coeffs, pixels, fill(order, length(pixels)), min_pixel, max_pixel, min_order, max_order)
    end
    return out
end

function polyfit2d(x, y, z, w, degx, degy)
    numpy = pyimport("numpy")
    V = numpy.polynomial.polynomial.polyvander2d(x, y, [degx, degy])
    Vw = V .* sqrt.(w[:, :])
    zw = z .* sqrt.(w)
    coeffs = Vw \ zw
    coeffs = collect(transpose(reshape(coeffs, degy + 1, degx + 1)))
    return coeffs
end

# Primary build method
function fit_peaks_poly2d(pixels, orders, λs, weights; deg_inter_order::Int, deg_intra_order::Int, nσ::Real=4)

    # Copy to not modify original
    pixels = copy(pixels)
    orders = copy(orders)
    λs = copy(λs)
    weights = copy(weights)

    # Actual vals to fit
    mλ = orders .* λs

    # Bounds for polynomials
    min_pixel = 1
    max_pixel = maximum(pixels)
    min_order = nanminimum(orders) - 1
    max_order = nanmaximum(orders) + 1

    # Get normalized vals
    pixels_norm = normalize_val.(pixels, min_pixel, max_pixel)
    orders_norm = normalize_val.(orders, min_order, max_order)

    # Good vals
    good = findall(isfinite.(weights))
    coeffs = nothing

    # Iteratively fit
    while true

        # Fit
        x = pixels_norm[good]
        y = orders_norm[good]
        z = mλ[good]
        w = weights[good]
        coeffs = polyfit2d(x, y, z, w, deg_intra_order, deg_inter_order)

        # Model and residuals
        model = build_λsol2d_poly_flat(coeffs, pixels, orders, min_pixel, max_pixel, min_order, max_order)
        residuals = λs .- model
        σ = robust_stddev(residuals, weights; nσ=4)
        good_new = findall(weights .> 0 .&& isfinite.(weights) .&& abs.(residuals) .< nσ * σ)
        if length(good) == length(good_new)
            break
        else
            good = good_new
        end
    end

    # Outputs
    model = build_λsol2d_poly_flat(coeffs, pixels, orders, min_pixel, max_pixel, min_order, max_order)
    residuals = λs .- model
    residuals .= δλ2δv.(residuals, λs)

    # Return
    return (;coeffs, good, model, residuals, min_pixel, max_pixel, min_order, max_order)

end


function normalize_val(x::Real, lo::Real, hi::Real)
    μ = (lo + hi) / 2
    s = hi - lo
    return 2 * (x - μ) / s
end

function weighted_quantile(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}; q::Real=0.5)
    good = findall(@. isfinite(x) && (w > 0) && isfinite(w))
    if length(good) > 0
        xx = @view x[good]
        ww = @view w[good]
        return quantile(xx, Weights(ww), q)
    else
        return NaN
    end
end



function robust_stddev(x::AbstractArray{<:Real}, w::Union{AbstractArray{<:Real}, Nothing}=nothing; nσ::Real=4)
    if isnothing(w)
        w = ones(size(x))
    end
    med = weighted_quantile(x, w)
    adevs = abs.(med .- x)
    mad = weighted_quantile(adevs, w)
    good = findall(adevs .< 1.4826 * mad * nσ)
    if length(good) > 1
        return @views weighted_stddev(x[good], w[good])
    else
        return NaN
    end
end


function weighted_stddev(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}; μ::Union{Real, Nothing}=nothing)
    good = findall(@. isfinite(x) && (w > 0) && isfinite(w))
    if length(good) == 0
        return NaN
    end
    xx = @view x[good]
    ww = w[good]
    ww ./= sum(ww)
    if isnothing(μ)
        μ = nansum(xx .* ww) / nansum(ww)
    end
    dev = xx .- μ
    bias_estimator = 1.0 - sum(ww.^2)
    σ = sqrt(sum(dev .^2 .* ww) / bias_estimator)
    return σ
end