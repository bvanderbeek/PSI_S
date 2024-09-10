# GEODESIC FUNCTIONS (depricates unregistered Geodesics package)
# 1. Do everything in radians! Note that 'sind(x)' is > 2x  slower than 'sin(deg2rad(x))'...
#    Optional argument to convert degrees/radians
# 2. Use 'isapprox(a, b; atol = 0, rtol = Δ)' to compare angles to 90.0! Reasonable choice
#    for Δ may be 10*eps()? Default is sqrt(eps()).

macro cartesian(θ::Union{Vector, Symbol, Expr}, φ::Union{Vector, Symbol, Expr}, r::Union{Vector, Symbol, Expr}) 
    esc(:( @.($r*cos($θ)*cos($φ)), @.($r*cos($θ)*sin($φ)), @.($r*sin($θ)) )) 
end
macro spherical(x::Union{Vector, Symbol, Expr}, y::Union{Vector, Symbol, Expr}, z::Union{Vector, Symbol, Expr}) 
    esc(:( @.(asin($z/sqrt($x^2+$y^2+$z^2))), @.(atan($y,$x)), @.(sqrt($x^2+$y^2+$z^2)) )) 
end

function geo_to_cartesian(latitude, longitude, r = 1.0)
    sinλ, cosλ = sincos(longitude)
    sinθ, cosθ = sincos(latitude)
    x = r*cosθ*cosλ
    y = r*cosθ*sinλ
    z = r*sinθ

    return x, y, z
end
function geo_to_cartesian(latitude::Array, longitude::Array, r::Array)
    x, y, z = zeros(size(latitude)), zeros(size(longitude)), zeros(size(r))
    for i in eachindex(latitude)
        x[i], y[i], z[i] = geo_to_cartesian(latitude[i], longitude[i], r[i])
    end

    return x, y, z
end

function forward_geodesic(ϕ₁, λ₁, Δ, α; tf_degrees::Bool = true)
    # Convert input from degrees to radians
    if tf_degrees
        K = π/180.0
        ϕ₁ = K*ϕ₁
        λ₁ = K*λ₁
        Δ = K*Δ
        α = K*α
    end
    # Trig computations
    sinϕ₁, cosϕ₁ = sincos(ϕ₁)
    sinΔ, cosΔ = sincos(Δ)
    sinα, cosα = sincos(α)
    # Destination
    sinϕ₂ = sinϕ₁*cosΔ + cosϕ₁*sinΔ*cosα
    λ₂ = λ₁ + atan(sinα*sinΔ*cosϕ₁, cosΔ-sinϕ₁*sinϕ₂)
    ϕ₂ = asin(sinϕ₂)
    # Solution for reference point at poles
    # At ϕ₁ = ±90° λ₂ is non-unique but can be deduced from α
    # by considering the direction α is pointing when you rotate
    # the sphere-centered cartesian x-coordinate which points
    # towards 0°N, 0°S to the poles.
    #
    # Assumes that λ₁ = 0° when α is calculated via inverse_geodesic,
    # however, this to could also be corrected for.
    if isapprox(ϕ₁, -0.5*π) # ϕ₁ == -0.5*π
        λ₂ = α
    elseif isapprox(ϕ₁, 0.5*π) # ϕ₁ == 0.5*π
        λ₂ = π - α
        if λ₂ > π
            λ₂ = λ₂ - 2.0*π
        end
    end
    # Convert output from radians to degrees
    if tf_degrees
        ϕ₂ /= K
        λ₂ /= K
    end

    return ϕ₂, λ₂ 
end
function forward_geodesic(ϕ₁, λ₁, Δ::AbstractArray, α; tf_degrees::Bool = true)
    # Pre-allocate output
    ϕ = similar(Δ)
    λ = similar(Δ)
    # Loop to retrieve geographic coordinates
    for i in eachindex(Δ)
        ϕ[i], λ[i] = forward_geodesic(ϕ₁, λ₁, Δ[i], α; tf_degrees = tf_degrees)
    end

    return ϕ, λ
end

function inverse_geodesic(ϕ₁, λ₁, ϕ₂, λ₂; tf_degrees::Bool = true)
    # Convert input from degrees to radians
    if tf_degrees
        K = π/180.0
        ϕ₁ = K*ϕ₁
        λ₁ = K*λ₁
        ϕ₂ = K*ϕ₂
        λ₂ = K*λ₂
    end
    # Coordinate differences
    Δϕ = ϕ₂ - ϕ₁
    Δλ = λ₂ - λ₁
    # Assign null longitude if at polls for consistent results
    if isapprox(abs(ϕ₁), 0.5*π) # abs(ϕ₁) == 0.5*π
        λ₁ = 0.0
        Δλ = λ₂
    end
    # Trig computations
    sinϕ₁, cosϕ₁ = sincos(ϕ₁)
    sinϕ₂, cosϕ₂ = sincos(ϕ₂)
    sinΔλ, cosΔλ = sincos(Δλ)
    # Haversine formula for distance between two points on a sphere
    h = sqrt( (sin(0.5*Δϕ)^2) + cosϕ₁*cosϕ₂*(sin(0.5*Δλ)^2) )
    h = min(h, 1.0)
    Δ = 2*asin(h)
    # Bearing formula
    α = atan( sinΔλ*cosϕ₂, cosϕ₁*sinϕ₂ - sinϕ₁*cosϕ₂*cosΔλ )
    # Convert output from radians to degrees
    if tf_degrees
        Δ /= K
        α /= K
    end

    return Δ, α
end