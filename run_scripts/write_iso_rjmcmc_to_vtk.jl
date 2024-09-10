using JLD
using Statistics
using LinearAlgebra
using StaticArrays
using WriteVTK

# START INPUTS
rϵ, rη, rγ = (1.0, 0.0, 0.0) # Define Thomsen parameter ratios...not relevant for isotropic models
# STOP INPUTS

# Coordinate System Conversion Functions from the PSI Tomography Package

function ecef_vector(east_north_elv::NTuple{3, T}, latitude, longitude) where {T}
    # ECEF components for vector at (0°N, 0°E)
    w = (east_north_elv[3], east_north_elv[1], east_north_elv[2])
    # Rotate to geographic position
    c = π/180.0
    R = rotation_matrix((-c*latitude, c*longitude), (2, 3))
    sx = R[1,1]*w[1] + R[1,2]*w[2] + R[1,3]*w[3]
    sy = R[2,1]*w[1] + R[2,2]*w[2] + R[2,3]*w[3]
    sz = R[3,1]*w[1] + R[3,2]*w[2] + R[3,3]*w[3]

    return sx, sy, sz
end

function rotation_matrix(α, n)
    sinα, cosα = sincos(α)
    if n == 1
        return @SMatrix [1.0 0.0 0.0; 0.0 cosα -sinα; 0.0 sinα cosα]
    elseif n == 2
        return @SMatrix [cosα 0.0 sinα; 0.0 1.0 0.0; -sinα 0.0 cosα]
    elseif n == 3
        return @SMatrix [cosα -sinα 0.0; sinα cosα 0.0; 0.0 0.0 1.0]
    else
        error("Requested rotation axis out-of-bounds!")
    end
end

function rotation_matrix(α::Tuple, n::Tuple)
    R = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    for i in eachindex(n)
        R = rotation_matrix(α[i],n[i])*R
    end
    return R
end

name = ARGS[1]
wrk_dir = string(@__DIR__, "/output/", name, "/")
gridded_model = string(wrk_dir, "/$(name)_grid_model.jld")
vtk_model = string(wrk_dir, "/$(name)_grid_model")

# Permute the models such that x1 = longitude, x2 = latitude, and x3 = Radial
Model = load(gridded_model);
dlnVp = permutedims(Model["dlnVp_mean"], (2, 1, 3));
rstd_dlnVp = permutedims(Model["dlnVp_relstd"], (2, 1, 3));

# Extract coordinate vectors
latitude = Model["latitude"]; # Radians
longitude = Model["longitude"]; # Radians
radius = Model["radius"]; # km

# Global coordinate arrays
nx, ny, nz = (length(longitude), length(latitude), length(radius))
xg, yg, zg = (zeros((nx, ny, nz)), zeros((nx, ny, nz)), zeros((nx, ny, nz)));
for (k, x3) in enumerate(radius)
    for (j, x2) in enumerate(latitude)
        sinϕ, cosϕ = sincos(x2)
        for (i, x1) in enumerate(longitude)
            sinλ, cosλ = sincos(x1)
            xg[i,j,k] = radius[k]*cosϕ*cosλ
            yg[i,j,k] = radius[k]*cosϕ*sinλ
            zg[i,j,k] = radius[k]*sinϕ
        end
    end
end

# Write VTK file
vtk_grid(vtk_model, xg, yg, zg) do vtk
    vtk["dlnVp"] = dlnVp
    vtk["rstd_dlnVp"] = rstd_dlnVp
end
