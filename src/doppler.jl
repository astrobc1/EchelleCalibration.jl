export δλ2δv, δv2δλ, δx2δλ

δλ2δv(δλ::Real, λ::Real) = δλ / λ * SPEED_OF_LIGHT_MPS
δv2δλ(δv::Real, λ::Real) = λ * (exp(δv / SPEED_OF_LIGHT_MPS) - 1)

function δx2δλ(δx::Real, x::Real, λsol::AbstractArray{<:Real})
    nx = length(λsol)
    xarr = collect(1:nx)
    good = findall(isfinite.(λsol))
    spl = @views DataInterpolations.CubicSpline(λsol[good], xarr[good])
    δλ = DataInterpolations.derivative(spl, x) * δx
    return δλ
end