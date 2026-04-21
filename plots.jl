using Random
using Plots
using LaTeXStrings
using GLM
using DataFrames
using Statistics
using CSV
using StatsBase
using LsqFit
using LinearAlgebra

default(xlabelfontsize=18, ylabelfontsize=18, xtickfontsize=14, ytickfontsize=14,
        left_margin=5Plots.mm, right_margin=5Plots.mm)

# ==============================================================================
# Paths
# ==============================================================================

const BASE_PATH   = raw"C:\Users\magic\OneDrive\Desktop\Programas\Julia_programs\programas_julia\ising_model_2d\jdatasets"
const COMBINED    = joinpath(BASE_PATH, "results_with_0.001_metrop_wolff_combinned")
const CRIT_PATH   = joinpath(BASE_PATH, "results_for_crit_region_L_256_for_crit_exponent_jump_0.001")
const CORR_PATH   = joinpath(BASE_PATH, "results_for_correlation_time_and_corr_length")

const L_LIST  = [4, 8, 16, 32, 64, 128]
const N_LIST  = L_LIST .^ 2          # number of spins
const COLORS  = [:blue, :purple, :green, :red, :orange, :yellow]
const N_MEAS  = 1000                 # number of measurements per simulation

# ==============================================================================
# Data loading
# ==============================================================================

"""Load correlation-time CSV for a given L from `dir`."""
function load_tau(L, dir=COMBINED)
    df = CSV.read(joinpath(dir, "correlation_time_L_$(L)_ising_2d.csv"), DataFrame)
    return (H=vec(df[:,1]), M=vec(df[:,2]), H2=vec(df[:,3]), M2=vec(df[:,4]), M4=vec(df[:,6]))
end

"""Load simulation-data CSV for a given L from `dir`."""
function load_measures(L, dir=COMBINED)
    df = CSV.read(joinpath(dir, "simulation_data_$(L)_ising_2d.csv"), DataFrame)
    return df
end

"""
Compute statistical error for observable X using integrated autocorrelation time τ:
    δX = σ_X * sqrt((2τ + 1) / N_meas)
"""
error_from_tau(var_X, tau) = sqrt.(var_X) .* sqrt.((2 .* tau .+ 1) ./ N_MEAS)

"""
Extract mean, variance, and autocorrelation-corrected error for H and M from a measures DataFrame.
Returns a NamedTuple with fields H, varH, M, varM, H2, varH2, M2, varM2, M4, varM4, T.
"""
function extract_observables(df, tau)
    H    = vec(df[:,1]);  varH    = vec(df[:,2]).^2
    M    = vec(df[:,3]);  varM    = vec(df[:,4]).^2
    H2   = vec(df[:,5]);  varH2   = vec(df[:,6]).^2
    M2   = vec(df[:,7]);  varM2   = vec(df[:,8]).^2
    M4   = vec(df[:,10]); varM4   = vec(df[:,11]).^2
    T    = vec(df[:,9])

    εH  = error_from_tau(varH,  tau.H)
    εM  = error_from_tau(varM,  tau.M)
    εH2 = error_from_tau(varH2, tau.H2)
    εM2 = error_from_tau(varM2, tau.M2)
    εM4 = error_from_tau(varM4, tau.M4)

    return (; H, varH, M, varM, H2, varH2, M2, varM2, M4, varM4, T, εH, εM, εH2, εM2, εM4)
end

"""Load all L data at once; returns a Dict keyed by L."""
function load_all(dir=COMBINED)
    T_ref = vec(CSV.read(joinpath(dir, "simulation_data_4_ising_2d.csv"), DataFrame)[:,9])
    data  = Dict{Int,Any}()
    for L in L_LIST
        tau = load_tau(L, dir)
        df  = load_measures(L, dir)
        obs = extract_observables(df, tau)
        data[L] = (; tau, obs, T=obs.T)
    end
    return data, T_ref
end

# ==============================================================================
# Physic quantities
# ==============================================================================

per_spin(x, N)  = x ./ N

energy_per_spin(obs, N)         = per_spin(obs.H, N), per_spin(obs.εH, N)
magnetisation_per_spin(obs, N)  = per_spin(obs.M, N), per_spin(obs.εM, N)

function heat_capacity(obs, N)
    cv  = obs.T.^(-2) .* obs.varH ./ N
    δcv = obs.T.^(-2) .* (obs.εH .+ obs.εH2) ./ N
    return cv, δcv
end

function susceptibility(obs, N)
    χ  = (1 ./ obs.T) .* obs.varM ./ N
    δχ = (1 ./ obs.T) .* (obs.εM .+ obs.εM2) ./ N
    return χ, δχ
end

function binder_cumulant(obs)
    U  = 1 .- obs.M4 ./ (3 .* obs.M2.^2)
    δU = (1 ./(3 .* obs.M2.^2)) .* obs.εM4 .+
         (2 .* obs.M4 ./ (3 .* obs.M2.^3)) .* obs.εM2
    return U, δU
end

# ==============================================================================
# Plot helpers
# ==============================================================================

"""Scatter all L series onto an existing plot object `p`."""
function plot_all_L!(p, T, ys, δys=nothing; kwargs...)
    for (i, L) in enumerate(L_LIST)
        yerr = isnothing(δys) ? nothing : δys[i]
        scatter!(p, T, ys[i]; yerror=yerr, label="L = $L", color=COLORS[i], kwargs...)
    end
end

save(name) = savefig(joinpath(@__DIR__, name))

# ==============================================================================
# Main analysis Functions
# ==============================================================================

function plot_τ()
    data_corr = Dict{Int,Any}()
    T_array   = 3.00:-0.01:2.00

    # Correlation-length data (only L = 64, 128 available)
    for L in [64, 128]
        df = CSV.read(joinpath(BASE_PATH, "correlation_function_$(L)_ising_2d.csv"), DataFrame)
        data_corr[L] = (ξ=collect(df[:,1]), σξ=collect(df[:,2]))
    end

    for L in [64, 128]
        d = data_corr[L]
        scatter(T_array, d.ξ; yerror=d.σξ, label="L = $L",
                xlabel=L"T", ylabel=L"\xi")
        save("corr_len_vs_T_L_$(L).png")
    end
end


function plot_measures()
    data, _ = load_all()

    # Build series vectors
    us = []; δus = []; ms = []; δms = []
    cvs = []; δcvs = []; χs = []; δχs = []
    Us = []; δUs = []

    T = data[4].T   # same T for all L

    for L in L_LIST
        obs = data[L].obs
        N   = L^2
        u, δu   = energy_per_spin(obs, N)
        m, δm   = magnetisation_per_spin(obs, N)
        cv, δcv = heat_capacity(obs, N)
        χ, δχ   = susceptibility(obs, N)
        U, δU   = binder_cumulant(obs)

        push!(us, u);   push!(δus, δu)
        push!(ms, m);   push!(δms, δm)
        push!(cvs, cv); push!(δcvs, δcv)
        push!(χs, χ);   push!(δχs, δχ)
        push!(Us, U);   push!(δUs, δU)
    end

    # Binder cumulant near Tc
    p_U = scatter(xlabel=L"T", ylabel=L"U", xlim=[2.255, 2.275])
    plot_all_L!(p_U, T, Us, δUs)
    save("binder_cumulant.png")
end


function crit_exponents()
    Tc = 2.269

    # ── L = 256 near critical region ──────────────────────────────────────────
    tau256  = load_tau(256, CRIT_PATH)
    df256   = load_measures(256, CRIT_PATH)

    t_min = findfirst(df256.T .== 2.0)
    t_max = findfirst(df256.T .< 2.267)
    df256   = df256[t_max:t_min, :]

    # Trim tau arrays to the same window
    tau256_trim = (; H=tau256.H[t_max:t_min], M=tau256.M[t_max:t_min],
                     H2=tau256.H2[t_max:t_min], M2=tau256.M2[t_max:t_min],
                     M4=tau256.M4[t_max:t_min])

    obs256 = extract_observables(df256, tau256_trim)
    T256   = obs256.T
    N256   = 256^2

    u256, _   = energy_per_spin(obs256, N256)
    m256, _   = magnetisation_per_spin(obs256, N256)
    cv256, _  = heat_capacity(obs256, N256)
    χ256, _   = susceptibility(obs256, N256)

    t_list = abs.((T256 ./ Tc) .- 1)

    # Log-log linear fits
    function log_fit(log_t, log_q, label)
        valid = isfinite.(log_t) .& isfinite.(log_q)
        df_fit = DataFrame(x=log_t[valid], y=log_q[valid])
        model  = lm(@formula(y ~ x), df_fit)
        n, m   = coef(model)
        δn, δm = stderror(model)
        R2     = r2(model)
        println("$label → intercept=$n ± $δn, slope=$m ± $δm, R²=$R2")
        return model, df_fit
    end

    model_m,  df_m  = log_fit(log.(t_list), log.(m256),  "m")
    model_cv, df_cv = log_fit(log.(t_list), log.(cv256), "cv")
    model_χ,  df_χ  = log_fit(log.(t_list), log.(χ256),  "χ")

    # ── χ_max ~ L^(γ/ν) fit across L sizes ───────────────────────────────────
    data, T_ref = load_all()
    χ_max = Float64[]
    for L in L_LIST
        obs = data[L].obs; N = L^2
        χ, _ = susceptibility(obs, N)
        push!(χ_max, maximum(χ))
    end

    model_χL, _ = log_fit(log.(Float64.(L_LIST)), log.(χ_max), "χ_max vs L")
    _, m_χL     = coef(model_χL)
    _, δm_χL    = stderror(model_χL)
    println("ν estimate: $(1.75/m_χL) ± $((1.75/m_χL^2)*δm_χL)")

    # ── Plots ─────────────────────────────────────────────────────────────────
    function fit_plot(df_fit, model, xlabel_str, ylabel_str)
        p = plot(xlabel=xlabel_str, ylabel=ylabel_str)
        scatter!(p, df_fit.x, df_fit.y; label="Data", color=:blue)
        plot!(p, df_fit.x, predict(model); label="Fit", color=:red, linestyle=:dash, lw=2)
        return p
    end

    p1 = fit_plot(df_m,  model_m,  L"\log|t|", L"\log(m)")
    p2 = fit_plot(df_cv, model_cv, L"\log|t|", L"\log(c_v)")
    p3 = fit_plot(df_χ,  model_χ,  L"\log|t|", L"\log(\chi)")

    log_χmax_df = DataFrame(x=log.(Float64.(L_LIST)), y=log.(χ_max))
    p4 = fit_plot(log_χmax_df, model_χL, L"\log(L)", L"\log(\chi_{\max})")

    combined = plot(p1, p2, p3, p4; layout=(2,2), size=(800,600), left_margin=5Plots.mm)
    display(combined); save("crit_exponents.png")
    for (p, name) in zip([p1,p2,p3,p4], ["m","cv","chi","chi_nu"])
        savefig(p, joinpath(@__DIR__, "crit_exponents_$(name).png"))
    end
end


function finite_size()
    Tc    = 2.269
    data, _ = load_all()
    T     = data[4].T
    t     = (T ./ Tc) .- 1

    # Known critical exponents for 2D Ising
    γ_over_ν = 1.75   # γ/ν
    one_over_ν = 1.0  # 1/ν  (ν = 1)

    p = plot(xlabel=L"t\, L^{1/\nu}", ylabel=L"\chi\, L^{-\gamma/\nu}",
             bottom_margin=5Plots.mm)

    for (i, L) in enumerate(L_LIST)
        obs = data[L].obs; N = L^2
        χ, _ = susceptibility(obs, N)
        scatter!(p, t .* L^one_over_ν, χ .* L^(-γ_over_ν);
                 label="L = $L", color=COLORS[i], marker=true)
    end

    display(p)
    save("finite_size_scaling.png")
end


# ==============================================================================
# Entry point
# ==============================================================================

# Uncomment the function(s) you want to run:
# plot_τ()
# plot_measures()
# crit_exponents()
finite_size()