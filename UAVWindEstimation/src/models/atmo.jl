struct Atmosphere{T, I1, I2}
    ρ::T
    lat::Vector{T}
    lon::Vector{T}
    u_wind_table::I1
    v_wind_table::I2
end

function gfs_wind_atmo(date_str, p_level)
    filepath = joinpath(DATADIR, date_str, "gfs_atmo.csv")
    df = CSV.read(filepath, DataFrame, header=["d1", "d2", "field", "ident", "lon", "lat", "val"])
    mask = endswith.(df.ident, "$p_level mb")
    df = df[mask, :]

    lat = unique(df.lat)
    lon = unique(df.lon)
    u_wind = [ first(df[(df.lon.==lon_j) .& (df.lat.==lat_i) .& (df.field.=="UGRD"), :]).val for lat_i in lat, lon_j in lon]
    v_wind = [ first(df[(df.lon.==lon_j) .& (df.lat.==lat_i) .& (df.field.=="VGRD"), :]).val for lat_i in lat, lon_j in lon]

    ρ = 0.41

    lat = Vector{Float64}(lat)
    lon = Vector{Float64}(lon)

    u_itp = LinearInterpolation((lat, lon), u_wind)
    v_itp = LinearInterpolation((lat, lon), v_wind)
    return Atmosphere(ρ, lat, lon, u_itp, v_itp)
end

compute_atmo_wind(a::Atmosphere, lat, lon) = a.ρ, a.u_wind_table(lat, lon), a.v_wind_table(lat, lon)

