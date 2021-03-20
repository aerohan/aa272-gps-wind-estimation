#
# Truth
#
struct EphemerisErrors{T}
    r_errors::Vector{SVector{3, T}}
    v_errors::Vector{SVector{3, T}}
    B_errors::Vector{T}
end
function EphemerisErrors(N, σ_pos::T, σ_vel::T, σ_B::T) where T<:Real
    r_errors = Vector{SVector{3, T}}()
    v_errors = Vector{SVector{3, T}}()
    B_errors = Vector{T}()

    for svidx=1:N
        r_error = (@SVector randn(3))*σ_pos
        v_error = (@SVector randn(3))*σ_vel
        B_error = randn()*σ_B
        push!(r_errors, r_error)
        push!(v_errors, v_error)
        push!(B_errors, B_error)
    end

    return EphemerisErrors(r_errors, v_errors, B_errors)
end

compute_pseudorange(r_rcv, r_sat, b_u, B_sat, ϵ_err) = norm(r_sat - r_rcv) + b_u + B_sat + ϵ_err
function compute_pseudorange_rate(r_rcv, r_sat, v_rcv, v_sat, ϵ_err)
    v_rel = v_sat - v_rcv
    r_rel = r_sat - r_rcv
    r̂_rel = r_rel/norm(r_rel)
    v_rel'r̂_rel + ϵ_err
end

function compute_pseudoranges_truth!(ρ_result, ρ̇_result, r_rcv, v_rcv, 
                                    eph_data::GpsEphemeris, errors::EphemerisErrors, b_u, σ_ρ, σ_ρ̇, gps_time)
    ρ_result .= NaN
    ρ̇_result .= NaN

    for svidx=1:length(ρ_result)
        r_eph, v_eph, B_sat = compute_ephemeris(gps_time, eph_data.data[svidx])

        r_sat_true = r_eph + errors.r_errors[svidx]
        v_sat_true = v_eph + errors.v_errors[svidx]

        if satellite_visible(r_rcv, r_sat_true)
            ρ_result[svidx] = compute_pseudorange(r_rcv, r_sat_true, b_u, B_sat, randn()*σ_ρ)
            ρ̇_result[svidx] = compute_pseudorange_rate(r_rcv, r_sat_true, v_rcv, v_sat_true, randn()*σ_ρ̇)
        end
    end

end

function satellite_visible(r_rcv, r_sat)
    r_rcv2sat_enu = ENU(ECEF(r_sat), ECEF(r_rcv), wgs84)
    elevation = atan(r_rcv2sat_enu[3], √(r_rcv2sat_enu[1]^2 + r_rcv2sat_enu[2]^2))
    return elevation*180/π >= 10.0
end

#
# Nav
#
function compute_pseudoranges_nav!(ρ_result, ρ̇_result, r_rcv, v_rcv, 
                                    eph_data::GpsEphemeris, b_u, gps_time)
    ρ_result .= NaN
    ρ̇_result .= NaN

    for svidx=1:length(ρ_result)
        r_eph, v_eph, B_sat = compute_ephemeris(gps_time, eph_data.data[svidx])

        ρ_result[svidx] = compute_pseudorange(r_rcv, r_eph, b_u, B_sat, 0.0)
        ρ̇_result[svidx] = compute_pseudorange_rate(r_rcv, r_eph, v_rcv, v_eph, 0.0)
    end

end

