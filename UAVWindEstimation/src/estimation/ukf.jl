α(i, κ, L) = i == 0 ? κ/(L+κ) : 1/(2*(L+κ))

function ukf_step(y_meas, μ_x_last, P_x_last, u, Q, R, κ, dynamics_func, measurement_func)
    nx = length(μ_x_last)
    ny = length(y_meas)

    μ_z_last = vcat(μ_x_last, @SVector(zeros(nx)))
    Σ_zz = [P_x_last 0I; 0I Q]
    L1 = cholesky(Symmetric(Σ_zz)).L

    z_sigma = typeof(μ_z_last)[]
    L = 2*nx
    push!(z_sigma, μ_z_last)
    for i=1:L
        push!(z_sigma, μ_z_last + √(L+κ)*L1[:,i])
        push!(z_sigma, μ_z_last - √(L+κ)*L1[:,i])
    end

    x_prop_sigma = map(z_sigma) do z
        x_prev_sigma = z[1:nx]
        w_sigma = z[nx+1:end]
        x_prop = dynamics_func(x_prev_sigma, u, w_sigma)
        x_prop
    end

    α_weights = [α(i, κ, L) for i=0:2*L]
    μ_x_prop = sum(α_weights .* x_prop_sigma)
    P_x_prop = sum([α_i * (x_prop_i .- μ_x_prop) * (x_prop_i .- μ_x_prop)' 
                    for (α_i, x_prop_i) in zip(α_weights, x_prop_sigma)])

    μ_z_prop = vcat(μ_x_prop, @SVector(zeros(ny)))
    Σ_zz = [P_x_prop zeros(nx, ny); zeros(ny, nx) R]
    L2 = cholesky(Symmetric(Σ_zz)).L

    z_sigma = typeof(μ_z_prop)[]
    L = nx + ny
    push!(z_sigma, μ_z_prop)
    for i=1:L
        push!(z_sigma, μ_z_prop + √(L+κ)*L2[:,i])
        push!(z_sigma, μ_z_prop - √(L+κ)*L2[:,i])
    end

    xy_pred_sigma = map(z_sigma) do z
        x_prop_sigma = z[1:nx]
        n_sigma = z[nx+1:end]
        y_meas_pred = measurement_func(x_prop_sigma, n_sigma)
        x_prop_sigma, y_meas_pred
    end

    α_weights = [α(i, κ, L) for i=0:2*L]
    μ_y_pred = sum(α_weights .* [y_pred_i for (_, y_pred_i) in xy_pred_sigma])
    Σ_yy = sum([α_i * (y_pred_i .- μ_y_pred) * (y_pred_i .- μ_y_pred)' 
                for (α_i, (_, y_pred_i)) in zip(α_weights, xy_pred_sigma)])
    Σ_xy = sum([α_i * (x_prop_i .- μ_x_prop) * (y_pred_i .- μ_y_pred)' 
                for (α_i, (x_prop_i, y_pred_i)) in zip(α_weights, xy_pred_sigma)])

    K = Σ_xy * inv(Symmetric(Σ_yy))
    P̂ = Symmetric(P_x_prop - K * Σ_xy')
    x̂ = μ_x_prop + K * (y_meas - μ_y_pred)

    x̂, P̂
end
