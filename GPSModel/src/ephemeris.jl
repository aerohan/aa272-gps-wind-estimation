struct SvEphemeris{T}
    sqrtA::T
    DeltaN::T
    Eccentricity::T
    Toe::T
    M0::T
    omega::T
    Cus::T
    Cuc::T
    Crs::T
    Crc::T
    Cis::T
    Cic::T
    Omega0::T
    OmegaDot::T
    Io::T
    IDOT::T
    SVclockBias::T
end

struct GpsEphemeris{T}
    data::Vector{SvEphemeris{T}}
end
Base.length(g::GpsEphemeris) = length(g.data)

function broadcast_ephemeris(date_str, time)
    gr = pyimport("georinex")
    eph_path = joinpath(DATADIR, date_str, "broadcast_eph.20n")
    eph = gr.load(eph_path)

    times = unix2datetime.(eph.time.data.astype("datetime64[s]").astype("int"))

    sqrtA_data = eph.sqrtA.data
    DeltaN_data = eph.DeltaN.data
    Eccentricity_data = eph.Eccentricity.data
    Toe_data = eph.Toe.data
    M0_data = eph.M0.data
    omega_data = eph.omega.data
    Cus_data = eph.Cus.data
    Cuc_data = eph.Cuc.data
    Crs_data = eph.Crs.data
    Crc_data = eph.Crc.data
    Cis_data = eph.Cis.data
    Cic_data = eph.Cic.data
    Omega0_data = eph.Omega0.data
    OmegaDot_data = eph.OmegaDot.data
    Io_data = eph.Io.data
    IDOT_data = eph.IDOT.data
    SVclockBias_data = eph.SVclockBias.data

    function get_single_eph_value(eph_vec)
        eph_vec[findfirst((.~isnan.(eph_vec)) .& (times .> time))]
    end

    all_eph_data = Vector{SvEphemeris{Float64}}()
    for sv = 1:length(eph.sv)
        sveph = SvEphemeris(
            get_single_eph_value(sqrtA_data[:,sv]),
            get_single_eph_value(DeltaN_data[:,sv]),
            get_single_eph_value(Eccentricity_data[:,sv]),
            get_single_eph_value(Toe_data[:,sv]),
            get_single_eph_value(M0_data[:,sv]),
            get_single_eph_value(omega_data[:,sv]),
            get_single_eph_value(Cus_data[:,sv]),
            get_single_eph_value(Cuc_data[:,sv]),
            get_single_eph_value(Crs_data[:,sv]),
            get_single_eph_value(Crc_data[:,sv]),
            get_single_eph_value(Cis_data[:,sv]),
            get_single_eph_value(Cic_data[:,sv]),
            get_single_eph_value(Omega0_data[:,sv]),
            get_single_eph_value(OmegaDot_data[:,sv]),
            get_single_eph_value(Io_data[:,sv]),
            get_single_eph_value(IDOT_data[:,sv]),
            get_single_eph_value(SVclockBias_data[:,sv]),
        )
        push!(all_eph_data, sveph) 
    end

    return GpsEphemeris(all_eph_data)
end

function compute_ephemeris(t, ephem)
    μ = 3.986004418e14
    Ω_e = 7.292115e-5

    a = (ephem.sqrtA)^2
    n = √(μ/a^3) + ephem.DeltaN
    e = ephem.Eccentricity

    tk = t - ephem.Toe
    Mk = ephem.M0 + n * tk
    result = nlsolve(n_ary(Ek->(Mk - Ek + e*sin(Ek))), [Mk])
    Ek = result.zero[1]

    sinνk = √(1 - e^2)*sin(Ek) / (1 - e*cos(Ek))
    cosνk = (cos(Ek) - e)/(1 - e*cos(Ek))
    νk = atan(sinνk, cosνk)

    ϕk = νk + ephem.omega
    δϕk = ephem.Cus*sin(2ϕk) + ephem.Cuc*cos(2ϕk)
    uk = ϕk + δϕk

    δrk = ephem.Crs*sin(2ϕk) + ephem.Crc*cos(2ϕk)
    δik = ephem.Cis*sin(2ϕk) + ephem.Cic*cos(2ϕk)

    Ωk = ephem.Omega0 - Ω_e*t + ephem.OmegaDot*tk

    rk = a*(1 - e*cos(Ek)) + δrk
    ik = ephem.Io + ephem.IDOT*tk + δik

    xp = rk * cos(uk)
    yp = rk * sin(uk)

    p = a*(1 - e^2)
    vxp = - sin(uk)*√(μ/p)
    vyp = (e + cos(uk))*√(μ/p)

    r_ecef = @SVector [
                        xp*cos(Ωk) - yp*cos(ik)*sin(Ωk),
                        xp*sin(Ωk) + yp*cos(ik)*cos(Ωk),
                        yp*sin(ik)
                   ]

    v_ecef = @SVector [
                        vxp*cos(Ωk) - vyp*cos(ik)*sin(Ωk),
                        vxp*sin(Ωk) + vyp*cos(ik)*cos(Ωk),
                        vyp*sin(ik)
                   ]

    sv_bias = ephem.SVclockBias * 299792458.0

    r_ecef, v_ecef, sv_bias
end

