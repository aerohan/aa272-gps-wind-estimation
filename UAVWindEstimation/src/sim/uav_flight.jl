#
# Simulation Model
#
@dynamics UAVDubinsDynamics{T} begin
    @integrable begin
        lat::T
        lon::T
        alt::T
        vmag_air::T
        yaw::T
    end

    @direct begin
        v_wind_ne::SVector{2, T}
        v_ne::SVector{2, T}
    end
end

function DynamicsAndControl.initialize(::Type{UAVDubinsDynamics}, config)
    @unpack lat, lon, alt, vmag_air, yaw = config.initial

    atmosphere = config.atmosphere

    return (lat, lon, alt, vmag_air, yaw), (SA[0.0, 0.0], SA[0.0, 0.0]), (;atmosphere)
end

function DynamicsAndControl.dynamics!(this::UAVDubinsDynamics, ẋ, x, u, t)
    @unpack lat, lon, alt, vmag_air, yaw = x

    ρ, u_wind, v_wind = compute_atmo_wind(static(this).atmosphere, lat*180/π, lon*180/π)
    x.v_wind_ne = SA[v_wind, u_wind]
    v_ne = SA[vmag_air*cos(yaw), vmag_air*sin(yaw)] + x.v_wind_ne
    x.v_ne = v_ne

    ẋ.lon = v_ne[2]/(6378e3*cos(lat))
    ẋ.lat = v_ne[1]/6378e3
    ẋ.alt = 0.0
    ẋ.vmag_air = 0.0
    ẋ.yaw = u.turn_rate
end

function DynamicsAndControl.update!(this::UAVDubinsDynamics, ẋ, x, u, t)
    @unpack lat, lon, alt, vmag_air, yaw, v_wind_ne = x
    v_ne = SA[vmag_air*cos(yaw), vmag_air*sin(yaw)] + v_wind_ne

    log!(this, :state, t, (;lat, lon, alt, vmag_air, yaw, v_wind_ne, v_ne))

    #r_nav_error =  r_ecef - u.r_est
    #log!(this, :nav_error, t, (;r=r_nav_error))

    println("\tt=$t")

    return false
end

@sensor UAVGnssSensor{T, N} begin
    @outputs begin
        ρ::MVector{N, T}
        ρ̇::MVector{N, T}
        vmag_air::T
        yaw::T
    end
end

function DynamicsAndControl.initialize(::Type{UAVGnssSensor}, config)
    gps_eph = broadcast_ephemeris(config.date_str, DateTime(0.0))
    N = length(gps_eph)

    @unpack σ_pos, σ_vel, σ_B, t0_gps, σ_ρ, σ_ρ̇, b_u = config
    eph_errors = EphemerisErrors(N, σ_pos, σ_vel, σ_B)
    static_conf = (ephemeris=gps_eph, eph_errors, t0_gps, b_u, σ_ρ, σ_ρ̇)

    outputs_init = (@MVector(zeros(N)), @MVector(zeros(N)), 0.0, 0.0)

    return (), outputs_init, static_conf
end

function DynamicsAndControl.update!(this::UAVGnssSensor, y, _, x, t)
    @unpack lat, lon, alt, v_ne = x
    r_rcv = SVector(ECEF(LLA(lat*180/π, lon*180/π, alt), wgs84))
    v_rcv = vne_to_ecef(v_ne, r_rcv)
    @unpack ephemeris, eph_errors, t0_gps, σ_ρ, σ_ρ̇, b_u = static(this)
    compute_pseudoranges_truth!(y.ρ, y.ρ̇, r_rcv, v_rcv, ephemeris, eph_errors, b_u, σ_ρ, σ_ρ̇, t + t0_gps)

    y.vmag_air = x.vmag_air + randn()*0.1
    y.yaw = x.yaw + randn()*0.1

    @unpack ρ, ρ̇ = deepcopy(y)
    log!(this, :gnss, t, (;ρ, ρ̇))
end

@controller UAVFlightComputer{T} begin
    @state begin
        μ::SVector{8, T}
        P::SMatrix{8, 8, T}

        turn_rate_prev::T

        periodic::PeriodicReal{T}
        yaw_des::T
    end

    @outputs begin
        # nav outputs
        lat_est::T
        lon_est::T
        alt_est::T
        vmag_air_est::T
        yaw_est::T
        v_wind_est::SVector{2, T}
        b_u_est::T

        # control outputs
        turn_rate::T
    end
end

function DynamicsAndControl.initialize(::Type{UAVFlightComputer}, config)
    lla_initial = LLA(config.lat0, config.lon0, config.alt0)
    r_ecef_initial = SVector(ECEF(lla_initial, wgs84))
    @unpack P0, t0_gps = config
    periodic = PeriodicReal(250.0)
    state_init = (SA[config.lat0*π/180*1e5, config.lon0*π/180*1e5, config.alt0, 10.0, 20*π/180, 0.0, 0.0, 0.0], P0, 0.0, periodic, 0.0)
    outputs_init = (config.lat0*π/180, config.lon0*π/180, config.alt0, 10.0, 20*π/180, @SVector(zeros(2)), 0.0, 0.0)

    @unpack Q, R_ρ, R_ρ̇, R_vmag, R_yaw = config
    @unpack K_guidance = config
    gps_eph = broadcast_ephemeris(config.date_str, DateTime(0.0))
    N = length(gps_eph)
    ρ_data = zeros(N)
    ρ̇_data = zeros(N)
    static_conf = (;Q, R_ρ, R_ρ̇, R_vmag, R_yaw, ρ_data, ρ̇_data, gps_eph, t0_gps, K_guidance)

    return state_init, outputs_init, static_conf
end

function DynamicsAndControl.update!(this::UAVFlightComputer, out, nav_state, y_sense, t)
    visible_mask = @. ~isnan(y_sense.ρ)

    dynamics_func = (xs, u, w) -> begin
        dt = 10.0
        lon = xs[1]/1e5
        lat = xs[2]/1e5
        alt = xs[3]
        vmag_air = xs[4]
        yaw = xs[5]
        v_wind_ne = xs[6:7]
        b_u = xs[8]

        v_ne = SA[vmag_air*cos(yaw), vmag_air*sin(yaw)] + v_wind_ne

        x_new = SA[
            1e5*(lon + v_ne[2]/(6378e3*cos(lat)) * dt),
            1e5*(lat + v_ne[1]/6378e3 * dt),
            alt,
            vmag_air,
            yaw + u * dt,
            v_wind_ne[1],
            v_wind_ne[2],
            b_u
        ] .+ w
        x_new
    end

    meas_func = (xs, n) -> begin
        @unpack ρ_data, ρ̇_data, gps_eph, t0_gps = static(this)
        r_ecef = SVector(ECEF(LLA(xs[1:2]*180/π/1e5..., xs[3]), wgs84))
        vmag_air = xs[4]
        yaw = xs[5]
        v_wind_ne = xs[6:7]
        v_ne = SA[vmag_air*cos(yaw), vmag_air*sin(yaw)] + v_wind_ne
        v_ecef = vne_to_ecef(v_ne, r_ecef)
        compute_pseudoranges_nav!(ρ_data, ρ̇_data, r_ecef, v_ecef, 
            gps_eph, xs[8], t + t0_gps)
        [ρ_data[visible_mask]..., ρ̇_data[visible_mask]..., xs[4], xs[5]] .+ n
    end

    @unpack Q, R_ρ, R_ρ̇, R_vmag, R_yaw = static(this)
    @unpack μ, P, yaw_des = nav_state
    y_sense_ukf = [y_sense.ρ[visible_mask]..., y_sense.ρ̇[visible_mask]..., y_sense.vmag_air, y_sense.yaw]
    R = collect(Diagonal([ones(count(visible_mask))*R_ρ..., ones(count(visible_mask))*R_ρ̇..., R_vmag, R_yaw]))
    μ, P = ukf_step(y_sense_ukf, μ, P, nav_state.turn_rate_prev, Q, R, 1.0, dynamics_func, meas_func)
    @pack! nav_state = μ, P

    lat_est = μ[1]/1e5
    lon_est = μ[2]/1e5
    alt_est = μ[3]
    vmag_air_est = μ[4]
    yaw_est = μ[5]
    v_wind_est = μ[6:7]
    b_u_est = μ[8]

    turn_rate = nav_state.turn_rate_prev
    if due!(nav_state.periodic, t)
        println("due at t=$t")
        turn_rate = randn()*0.05
        nav_state.turn_rate_prev = turn_rate
    end

    @pack! out = lat_est, lon_est, alt_est, vmag_air_est, yaw_est, b_u_est, turn_rate

    log!(this, :nav, t, (;lat_est, lon_est, alt_est, vmag_air_est, yaw_est, b_u_est, v_wind_est))
    log!(this, :guidance, t, (prev=nav_state.turn_rate_prev, turn_rate))
end

#
# Helper functions
#
function vne_to_ecef(vne, r_ecef)
    transform = ECEFfromENU(ECEF(r_ecef), wgs84)
    n̂ = transform(ENU(SA[0, 1, 0])) - r_ecef
    ê = transform(ENU(SA[1, 0, 0])) - r_ecef
    hcat(n̂, ê)*vne
end

#
# Simulation Cases and Plots
#
function DynamicsAndControl.simulate(::Type{UAVDubinsDynamics}, ::Val{:flight1})
    stanford_coords = 37.42826804631273, -122.17010639325883
    dynamics_init = (lat=stanford_coords[1]*π/180, lon=stanford_coords[2]*π/180, alt=10e3, vmag_air=12.0, yaw=30*π/180)
    dyn_config = (initial=dynamics_init, atmosphere=gfs_wind_atmo("2020-02-10", 200))

    b_u = 12.0
    t0_gps = 12312412.0
    date_str="2020-02-10"

    sense_config = (σ_pos=5.0, σ_vel=1.0, σ_B=3.0, σ_ρ=5.0, σ_ρ̇=0.1, t0_gps, b_u, date_str)
    nav_init = (lat0=stanford_coords[1]+.1, lon0=stanford_coords[2]-.2, alt0=10e3-121.0)
    P0 = SMatrix{8,8}(collect(Diagonal([.2*π/180*1e5, .2*π/180*1e5, 20, .1, 10π/180, 20.0, 20.0, 5.0])))
    Q = collect(Diagonal([100, 100, 1e-2, 0.01, 0.001, 0.1, 0.1, 0.1]))
    R_ρ = 10.0
    R_ρ̇ = 1.0
    R_vmag = 1.0
    R_yaw = 1.0
    K_guidance = 1.0
    nav_config = (;nav_init..., K_guidance, P0, date_str, t0_gps, Q, R_ρ, R_ρ̇, R_vmag, R_yaw)

    sim = Simulation(
        ( :truth, UAVDubinsDynamics, dyn_config ),
        ( :sensor, UAVGnssSensor, sense_config ),
        ( :control, UAVFlightComputer, nav_config ),
        3*3600.0, RK4(), dt=10.0
    )
    data = simulate(sim)
    return data
end

save_figs = false
const FIGPATH = joinpath(@__DIR__, "..", "..", "..", "figs")
function finalize(fname)
    save_figs && savefig(joinpath(FIGPATH, "lla.png"))
    plot!()
end

function RecipesBase.plot(data, ::Type{UAVDubinsDynamics}, ::Val{:lla})
    s = data.truth.state
    plot(
        plot(s.lat, s.lon, legend=false),
        plot(s.time, s.alt, legend=false),
        plot(s.time, s.v_wind_ne, legend=false),
        plot(s.time, s.yaw.*180.0/π, legend=false),
        layout=(2,2), size=(1200, 1200)
    )
    finalize("lla.png")
end

function RecipesBase.plot(data, ::Type{UAVDubinsDynamics}, ::Val{:prange})
    s = data.sensor.gnss
    plot(s.time, s.ρ[1])
    plot!(s.time, s.ρ[2])
    plot!(s.time, s.ρ[3])
end

function RecipesBase.plot(data, ::Type{UAVDubinsDynamics}, ::Val{:nav})
    s = data.truth.state
    n = data.control.nav
    #ne = data.truth.nav_error
    plot(
        begin
            plot(s.time, s.lat, label="truth")
            plot!(n.time, n.lat_est, label="lat est")
            plot!(ylabel="lat [deg]")
        end,
        begin
            plot(s.time, s.lon, label="truth")
            plot!(n.time, n.lon_est, label="lon est")
            plot!(ylabel="lon [deg]")
        end,
        begin
            plot(s.time, s.alt, label="truth")
            plot!(n.time, n.alt_est, label="alt est")
            plot!(ylabel="alt [m]")
        end,
        begin
            plot(s.time, s.v_wind_ne[1], label="n truth")
            plot!(n.time, n.v_wind_est[1], label="n est")
            plot!(s.time, s.v_wind_ne[2], label="e truth")
            plot!(n.time, n.v_wind_est[2], label="e est")
            plot!(ylabel="v wind [m/s]")
        end,
        begin
            plot(s.time, s.yaw, label="truth")
            plot!(n.time, n.yaw_est, label="est")
            plot!(ylabel="yaw angle [rad]")
        end,
        begin
            plot(s.time, s.vmag_air, label="vmag air")
            plot!(n.time, n.vmag_air_est, label="vmag air est")
            plot!(ylabel="airspeed [m/s]")
        end,
        layout=(3, 2), size=(1200, 1600), legend=(0.8, 0.6),
        xlabel="time [s]", margin=10mm
    )
end

function RecipesBase.plot(data, ::Type{UAVDubinsDynamics}, ::Val{:windmap})
    atmo = gfs_wind_atmo("2020-02-10", 200)
    x = [lon for lat in atmo.lat for lon in atmo.lon]
    y = [lat for lat in atmo.lat for lon in atmo.lon]
    scale = 0.03
    itp(lon, lat) = [atmo.u_wind_table(lat, lon), atmo.v_wind_table(lat, lon)].*scale
    quiver(x, y, quiver=itp, size=(800, 800))

    let s=data.truth.state
        plot!(s.lon.*180/π, s.lat.*180/π, color=:green, lw=10, label="trajectory")
        plot!(xlabel="longitude [deg]", ylabel="latitude [deg]")
    end
    finalize("windmap.png")
end

function RecipesBase.plot(data, ::Type{UAVDubinsDynamics}, ::Val{:bay})
    coast(region=[-125 -115 32 42], proj=:Mercator, land=:lightbrown, ocean=:seashell, figsize=12, show=true)
end

function RecipesBase.plot(data, ::Type{UAVDubinsDynamics}, ::Val{:guidance})
    s = data.truth.state
    g = data.control.guidance

    plot(g.time, g.turn_rate, label="turn_rate")
    plot!(g.time, g.prev, label="prev")
end


export UAVDubinsDynamics
