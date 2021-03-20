#
# Simulation Model
#
@dynamics StaticGnssDynamics{T} begin
    @integrable begin
        r_ecef::SVector{3,T}
    end
end

function DynamicsAndControl.initialize(::Type{StaticGnssDynamics}, config)
    lla_initial = LLA(config.lat0, config.lon0, config.alt0)
    r_ecef_initial = ECEF(lla_initial, wgs84)
    return (SVector(r_ecef_initial),), (), (;)
end

function DynamicsAndControl.dynamics!(this::StaticGnssDynamics, ẋ, x, u, t)
    ẋ.r_ecef = @SVector zeros(3)
end

function DynamicsAndControl.update!(this::StaticGnssDynamics, ẋ, x, u, t)
    @unpack r_ecef = x
    @unpack lat, lon, alt = LLA(ECEF(r_ecef), wgs84)

    log!(this, :state, t, (;r_ecef, lat, lon, alt))

    r_nav_error =  r_ecef - u.r_est
    log!(this, :nav_error, t, (;r=r_nav_error))

    println("\tt=$t")

    return false
end

@sensor StaticGnssSensor{T, N} begin
    @outputs begin
        ρ::MVector{T, N}
        ρ̇::MVector{T, N}
    end
end

function DynamicsAndControl.initialize(::Type{StaticGnssSensor}, config)
    gps_eph = broadcast_ephemeris(config.date_str, DateTime(0.0))
    N = length(gps_eph)

    @unpack σ_pos, σ_vel, σ_B, t0_gps, σ_ρ, b_u = config
    eph_errors = EphemerisErrors(N, σ_pos, σ_vel, σ_B)
    static_conf = (ephemeris=gps_eph, eph_errors, t0_gps, b_u, σ_ρ)

    outputs_init = (@MVector(zeros(N)), @MVector(zeros(N)))

    return (), outputs_init, static_conf
end

function DynamicsAndControl.update!(this::StaticGnssSensor, y, _, x, t)
    r_rcv = x.r_ecef
    v_rcv = @SVector zeros(3)
    @unpack ephemeris, eph_errors, t0_gps, σ_ρ, b_u = static(this)
    compute_pseudoranges_truth!(y.ρ, y.ρ̇, r_rcv, v_rcv, ephemeris, eph_errors, b_u, σ_ρ, t + t0_gps)

    @unpack ρ, ρ̇ = deepcopy(y)
    log!(this, :gnss, t, (;ρ, ρ̇))
end

@controller StaticNavComputer{T} begin
    @state begin
        μ::SVector{4, T}
        P::SMatrix{4, 4, T}
    end

    @outputs begin
        r_est::SVector{3, T}
        b_u_est::T
    end
end

function DynamicsAndControl.initialize(::Type{StaticNavComputer}, config)
    lla_initial = LLA(config.lat0, config.lon0, config.alt0)
    r_ecef_initial = SVector(ECEF(lla_initial, wgs84))
    @unpack P0, t0_gps = config
    state_init = (vcat(r_ecef_initial, @SVector([0.0])), P0)

    outputs_init = (r_ecef_initial, 0.0)

    #@unpack Q, R = config
    gps_eph = broadcast_ephemeris(config.date_str, DateTime(0.0))
    N = length(gps_eph)
    ρ_data = zeros(N)
    ρ̇_data = zeros(N)
    Q = Diagonal([10.0, 10.0, 10.0, 0.1])
    R = 1e5I
    static_conf = (;Q, R, ρ_data, ρ̇_data, gps_eph, t0_gps)

    return state_init, outputs_init, static_conf
end

function DynamicsAndControl.update!(this::StaticNavComputer, nav_out, nav_state, y_sense, t)
    visible_mask = @. ~isnan(y_sense.ρ)

    dynamics_func = (xs, u, w) -> xs .+ w
    meas_func = (xs, n) -> begin
        @unpack ρ_data, ρ̇_data, gps_eph, t0_gps = static(this)
        compute_pseudoranges_nav!(ρ_data, ρ̇_data, xs[1:3], @SVector(zeros(3)), 
            gps_eph, xs[4], t + t0_gps)
        ρ_data[visible_mask] .+ n
    end

    @unpack Q, R = static(this)
    @unpack μ, P = nav_state
    y_sense_pruned = y_sense.ρ[visible_mask]
    μ, P = ukf_step(y_sense_pruned, μ, P, nothing, Q, R, 1.0, dynamics_func, meas_func)
    @pack! nav_state = μ, P

    r_est = μ[1:3]
    b_u_est = μ[4]

    @pack! nav_out = r_est, b_u_est

    log!(this, :est, t, (;r_est, b_u_est))
end

#
# Simulation Cases and Plots
#
function DynamicsAndControl.simulate(::Type{StaticGnssDynamics}, ::Val{:simple_static})
    stanford_coords = 37.42826804631273, -122.17010639325883
    dyn_config = (lat0=stanford_coords[1], lon0=stanford_coords[2], alt0=5e3)

    b_u = 12.0
    t0_gps = 12312412.0
    date_str="2020-02-10"

    sense_config = (σ_pos=15.0, σ_vel=5.0, σ_B=7.0, σ_ρ=2.0, t0_gps, b_u, date_str)
    nav_init = (lat0=stanford_coords[1]+.1, lon0=stanford_coords[2]-.2, alt0=3e3)
    nav_config = (nav_init..., P0=SMatrix{4,4}(1e3I), date_str, t0_gps)

    sim = Simulation(
        ( :truth, StaticGnssDynamics, dyn_config ),
        ( :sensor, StaticGnssSensor, sense_config ),
        ( :nav, StaticNavComputer, nav_config ),
        360.0, RK4(), dt=1.0
    )
    data = simulate(sim)
    return data
end

function RecipesBase.plot(data, ::Type{StaticGnssDynamics}, ::Val{:lla})
    s = data.truth.state
    plot(
        plot(s.lat, s.lon, legend=false),
        plot(s.time, s.alt, legend=false),
        layout=(1,2), size=(1200, 600)
    )
end

function RecipesBase.plot(data, ::Type{StaticGnssDynamics}, ::Val{:prange})
    s = data.sensor.gnss
    plot(s.time, s.ρ[1])
    plot!(s.time, s.ρ[2])
    plot!(s.time, s.ρ[3])
end

function RecipesBase.plot(data, ::Type{StaticGnssDynamics}, ::Val{:nav})
    n = data.nav.est
    ne = data.truth.nav_error
    plot(
        plot(ne.time, ne.r),
        plot(ne.time, norm.(ne.r)),
        plot(n.time, n.b_u_est),
        layout=(1,3), size=(1200, 600)
    )
end

export StaticGnssDynamics
