#!/usr/bin/env julia

# ============================================================
# koala_wiser_null.jl
# ------------------------------------------------------------
# 作用：
#   1) 读取 R 端写出的 long-format 纵向数据；
#   2) 用 TrajGWAS/WiSER 拟合不含遗传效应的 null model；
#   3) 按 KOALA 推导，计算 mean 与 WS variability 两条 channel
#      所需的 subject-level (c, D, K, M)；
#   4) 将结果写回磁盘，供 R 端构造 augmented weighted LASSO。
#
# 输入 CSV 的列约定：
#   id, y,
#   X_1, ..., X_qx,          # mean fixed-effect covariates（不含截距）
#   Z_1, ..., Z_qz,          # random-effect covariates（不含截距）
#   W_1, ..., W_qw           # WSV fixed-effect covariates（不含截距）
#
# 注意：三条公式均由 Julia 自动加入截距 1。
# ============================================================

using CSV
using DataFrames
using LinearAlgebra
using DelimitedFiles
using TrajGWAS
using WiSER
using NLopt
using Serialization

BLAS.set_num_threads(1)

# -----------------------------
# SPA 扫描模式
# -----------------------------
# null-model 模式会将 TrajGWAS 的 invariant score-test object 与经验 CGF
# 保存为一个 Julia serialization 文件。R 端随后把同一 block 的 original
# 与 knockoff 矩阵一起送回本脚本；二者共用完全相同的 null model、ECGF、
# SPA cutoff 和数值流程。
function read_float64_matrix(path::String, n::Int, p::Int)
    isfile(path) || error("找不到 binary genotype 文件: $(path)")
    n > 0 || error("binary genotype 的 n 必须为正整数。")
    p > 0 || error("binary genotype 的 p 必须为正整数。")

    expected_bytes = Int64(n) * Int64(p) * Int64(sizeof(Float64))
    observed_bytes = filesize(path)
    observed_bytes == expected_bytes || error(
        "binary genotype 文件大小错误: file=$(path), " *
        "expected=$(expected_bytes), observed=$(observed_bytes)"
    )

    X = Matrix{Float64}(undef, n, p)
    open(path, "r") do io
        read!(io, X)
    end
    return X
end

function run_spa_mode()
    length(ARGS) >= 7 || error(
        "用法: julia koala_wiser_null.jl --spa-binary " *
        "<spa_state.jls> <G_original.f64> <G_knockoff.f64> " *
        "<output.csv> <n> <p> [r]"
    )

    state_file = ARGS[2]
    original_file = ARGS[3]
    knockoff_file = ARGS[4]
    output_file = ARGS[5]
    n = parse(Int, ARGS[6])
    p = parse(Int, ARGS[7])
    spa_r = length(ARGS) >= 8 ? parse(Float64, ARGS[8]) : 2.0

    isfile(state_file) || error("找不到 TrajGWAS SPA state: $(state_file)")
    isfile(original_file) || error("找不到 original binary genotype: $(original_file)")
    isfile(knockoff_file) || error("找不到 knockoff binary genotype: $(knockoff_file)")
    isfinite(spa_r) && spa_r >= 0 || error("SPA cutoff r 必须为有限非负数。")

    state = deserialize(state_file)
    st = state.scoretest
    Ks = state.ecgf
    expected_ids = String.(state.subject_ids)

    length(expected_ids) == n || error(
        "binary genotype 样本数与 SPA state 不一致: " *
        "binary n=$(n), state n=$(length(expected_ids))"
    )
    G = read_float64_matrix(original_file, n, p)
    G_tilde = read_float64_matrix(knockoff_file, n, p)

    P_beta = ones(Float64, p)
    P_tau = ones(Float64, p)
    dir_beta = zeros(Int, p)
    dir_tau = zeros(Int, p)
    P_chisq_beta = ones(Float64, p)
    P_chisq_tau = ones(Float64, p)

    P_tilde_beta = ones(Float64, p)
    P_tilde_tau = ones(Float64, p)
    dir_tilde_beta = zeros(Int, p)
    dir_tilde_tau = zeros(Int, p)
    P_tilde_chisq_beta = ones(Float64, p)
    P_tilde_chisq_tau = ones(Float64, p)

    function scan_one(g::AbstractVector{Float64})
        all(isfinite, g) || error("SPA genotype 中存在 NA/NaN/Inf。")

        # 与 TrajGWAS genome scan 对单态变异的处理保持一致。
        if maximum(g) == minimum(g)
            return 1.0, 1.0, 0, 0, 1.0, 1.0
        end

        p_beta_chisq, p_tau_chisq, d_beta_chisq, d_tau_chisq =
            test!(st, g, g)

        p_beta_spa, p_tau_spa, d_beta_spa, d_tau_spa = TrajGWAS.spa(
            g,
            st,
            (p_beta_chisq, p_tau_chisq),
            (d_beta_chisq, d_tau_chisq),
            Ks;
            r = spa_r
        )

        return p_beta_spa, p_tau_spa, d_beta_spa, d_tau_spa,
               p_beta_chisq, p_tau_chisq
    end

    for j in 1:p
        # TrajGWAS test! has historically expected a concrete Vector.
        # Copy one column at a time to preserve that API without materializing
        # an additional full genotype matrix.
        g = Vector{Float64}(@view G[:, j])
        g_tilde = Vector{Float64}(@view G_tilde[:, j])

        P_beta[j], P_tau[j], dir_beta[j], dir_tau[j],
        P_chisq_beta[j], P_chisq_tau[j] = scan_one(g)

        P_tilde_beta[j], P_tilde_tau[j],
        dir_tilde_beta[j], dir_tilde_tau[j],
        P_tilde_chisq_beta[j], P_tilde_chisq_tau[j] = scan_one(g_tilde)
    end

    result = DataFrame(
        index = collect(1:p),
        P_beta = P_beta,
        P_tau = P_tau,
        dir_beta = dir_beta,
        dir_tau = dir_tau,
        P_chisq_beta = P_chisq_beta,
        P_chisq_tau = P_chisq_tau,
        P_tilde_beta = P_tilde_beta,
        P_tilde_tau = P_tilde_tau,
        dir_tilde_beta = dir_tilde_beta,
        dir_tilde_tau = dir_tilde_tau,
        P_tilde_chisq_beta = P_tilde_chisq_beta,
        P_tilde_chisq_tau = P_tilde_chisq_tau
    )
    CSV.write(output_file, result)
    println("TrajGWAS SPA p values written to: ", output_file)
end

if !isempty(ARGS) && ARGS[1] == "--spa-binary"
    run_spa_mode()
    exit()
end

# -----------------------------
# 1. 读取命令行参数
# -----------------------------
if length(ARGS) < 2
    error("用法: julia koala_wiser_null.jl <input.csv> <output_dir> [runs] [maxeval]")
end

input_file = ARGS[1]
out_dir    = ARGS[2]
runs       = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10
maxeval    = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 4000

mkpath(out_dir)

# -----------------------------
# 2. 工具函数
# -----------------------------

"""按数字后缀排序列名，例如 X_1, X_2, ..., X_10。"""
function numbered_columns(df::DataFrame, prefix::String)
    cols = Symbol[]
    for nm in names(df)
        s = String(nm)
        if startswith(s, prefix)
            push!(cols, Symbol(s))
        end
    end
    sort!(cols, by = x -> parse(Int, split(String(x), "_")[end]))
    return cols
end

"""由列名动态构造 @formula(y ~ 1 + ...)。"""
function make_formula(cols::Vector{Symbol})
    rhs = isempty(cols) ? "1" : "1 + " * join(string.(cols), " + ")
    return eval(Meta.parse("@formula(y ~ $(rhs))"))
end

"""
构造 Σγ 的 vech 参数对应的对称基矩阵。
顺序为：(1,1), (2,1), ..., (q,1), (2,2), ..., (q,q)。
对于非对角元素，一个参数同时控制 Σ[a,b] 与 Σ[b,a]。
"""
function covariance_basis(q::Int)
    basis = Matrix{Float64}[]
    labels = String[]
    for b in 1:q
        for a in b:q
            B = zeros(Float64, q, q)
            if a == b
                B[a, b] = 1.0
            else
                B[a, b] = 1.0
                B[b, a] = 1.0
            end
            push!(basis, B)
            push!(labels, "Sigma_$(a)_$(b)")
        end
    end
    return basis, labels
end

"""写矩阵/向量到无表头 TSV。"""
function write_plain(path::String, x)
    open(path, "w") do io
        writedlm(io, x, '\t')
    end
end

# -----------------------------
# 3. 读入 long-format 数据并构造公式
# -----------------------------
df = CSV.read(input_file, DataFrame)

required = [:id, :y]
for nm in required
    nm in propertynames(df) || error("输入文件缺少列：$(nm)")
end

xcols = numbered_columns(df, "X_")
zcols = numbered_columns(df, "Z_")
wcols = numbered_columns(df, "W_")

mean_formula = make_formula(xcols)
re_formula   = make_formula(zcols)
wsv_formula  = make_formula(wcols)

println("Mean formula: ", mean_formula)
println("Random-effect formula: ", re_formula)
println("WSV formula: ", wsv_formula)

# -----------------------------
# 4. 拟合 WiSER / TrajGWAS null model
# -----------------------------
# 这里沿用用户参考脚本中的 NLopt + LD_MMA 设置。
solver = NLopt.Optimizer()
solver_config = Dict(
    "algorithm" => :LD_MMA,
    "maxeval"   => maxeval
)

nm = trajgwas(
    mean_formula,
    re_formula,
    wsv_formula,
    :id,
    df,
    nothing;
    solver = solver,
    solver_config = solver_config,
    runs = runs,
    nullfile = joinpath(out_dir, "wiser_null.txt")
)

println(nm)

# -----------------------------
# 5. 提取 null estimates
# -----------------------------
beta_hat  = Vector{Float64}(nm.β)
tau_hat   = Vector{Float64}(nm.τ)
Sigma_hat = Matrix{Float64}(nm.Σγ)

m = length(nm.data)
q_beta = length(beta_hat)           # mean nuisance dimension（含截距）
q_w    = length(tau_hat)            # WSV fixed nuisance dimension（含截距）
q_z    = size(Sigma_hat, 1)         # random-effect design dimension（含截距）

Sigma_basis, Sigma_labels = covariance_basis(q_z)
q_sigma = length(Sigma_basis)
r_tau   = q_w + q_sigma             # η = (τ, vech(Σγ))

# -----------------------------
# 6. 初始化 mean channel 的 c, D, K, M
# -----------------------------
c_beta = zeros(Float64, m)
d_beta = zeros(Float64, m)
K_beta = zeros(Float64, m, q_beta)
M_beta = zeros(Float64, q_beta, q_beta)

# -----------------------------
# 7. 初始化 WS variability channel 的 c, D, K, M
# -----------------------------
c_tau = zeros(Float64, m)
d_tau = zeros(Float64, m)
K_tau = zeros(Float64, m, r_tau)
M_tau = zeros(Float64, r_tau, r_tau)

# -----------------------------
# 8. 按受试者计算局部 score / expected curvature quantities
# -----------------------------
for i in 1:m
    obs = nm.data[i]

    # WiSER 内部已经按 subject 保存每个纵向块的设计矩阵。
    Xi = Matrix{Float64}(transpose(obs.Xt))
    Zi = Matrix{Float64}(transpose(obs.Zt))
    Wi = Matrix{Float64}(transpose(obs.Wt))
    yi = Vector{Float64}(obs.y)

    ni = length(yi)
    one_i = ones(Float64, ni)

    # ----- Null covariance -----
    # E_i = diag{exp(W_i τ_hat)}
    e_i = exp.(Wi * tau_hat)
    E_i = Diagonal(e_i)

    # V_i = E_i + Z_i Σγ Z_i^T
    V_i = Matrix(E_i) + Zi * Sigma_hat * transpose(Zi)
    V_i = 0.5 .* (V_i + transpose(V_i))

    # A_i = V_i^{-1}; 对每个 subject 的 ni 通常较小，直接通过 Cholesky 求逆。
    F_i = cholesky(Symmetric(V_i); check = true)
    A_i = Matrix(inv(F_i))
    A_i = 0.5 .* (A_i + transpose(A_i))

    # null residual r_bi = y_i - X_i beta_hat
    r_i  = yi - Xi * beta_hat
    R_bi = r_i * transpose(r_i) - V_i

    # ========================================================
    # Mean channel
    # c_beta,i = 1_i^T A_i r_i
    # d_beta,i = 1_i^T A_i 1_i
    # K_beta,i = 1_i^T A_i X_i
    # M_beta   = Σ X_i^T A_i X_i
    # ========================================================
    c_beta[i] = dot(one_i, A_i * r_i)
    d_beta[i] = dot(one_i, A_i * one_i)
    K_beta[i, :] .= vec(transpose(one_i) * A_i * Xi)
    M_beta .+= transpose(Xi) * A_i * Xi

    # ========================================================
    # WS variability channel
    # c_tau,i = tr(A_i R_bi A_i E_i)
    # d_tau,i = tr(A_i E_i A_i E_i)
    # ========================================================
    c_tau[i] = tr(A_i * R_bi * A_i * E_i)
    d_tau[i] = tr(A_i * E_i * A_i * E_i)

    # 构造 η=(τ, vech(Σγ)) 的全部 V_{i,a} = ∂V_i/∂η_a
    Vdirs = Vector{Matrix{Float64}}(undef, r_tau)

    # 对 τ_a：V_{i,a} = diag(e_i .* W_i[:,a])
    for a in 1:q_w
        Vdirs[a] = Matrix(Diagonal(e_i .* Wi[:, a]))
    end

    # 对 s_r = vech(Σγ)：V_{i,r} = Z_i L_r Z_i^T
    for r in 1:q_sigma
        Vdirs[q_w + r] = Zi * Sigma_basis[r] * transpose(Zi)
    end

    # K_tau[i,a] = tr(A_i E_i A_i V_{i,a})
    for a in 1:r_tau
        K_tau[i, a] = tr(A_i * E_i * A_i * Vdirs[a])
    end

    # M_tau[a,b] = Σ_i tr(A_i V_{i,a} A_i V_{i,b})
    # 这里只累加上三角，再镜像到下三角。
    AV = [A_i * Vdirs[a] for a in 1:r_tau]
    for a in 1:r_tau
        for b in a:r_tau
            val = tr(AV[a] * AV[b])
            M_tau[a, b] += val
            if a != b
                M_tau[b, a] += val
            end
        end
    end
end

# 数值对称化
M_beta = 0.5 .* (M_beta + transpose(M_beta))
M_tau  = 0.5 .* (M_tau  + transpose(M_tau))

# -----------------------------
# 9. 输出结果
# -----------------------------
# nm.ids 在参考 TrajGWAS 代码中用于恢复 subject 顺序。
subject_ids = unique(string.(nm.ids))
length(subject_ids) == m || error("nm.ids 与 nm.data 的 subject 数量不一致。")

# TrajGWAS 原生 empirical SPA state。WSVarScoreTestInvariant(nm, 1, 1)
# 与 trajgwas(nm, geno_path; usespa=true) 的 single-SNP 路径完全一致。
spa_scoretest = WSVarScoreTestInvariant(nm, 1, 1)
spa_ecgf = TrajGWAS.ecgf(spa_scoretest)
serialize(
    joinpath(out_dir, "trajgwas_spa_state.jls"),
    (
        scoretest = spa_scoretest,
        ecgf = spa_ecgf,
        subject_ids = subject_ids
    )
)

write_plain(joinpath(out_dir, "subject_ids.tsv"), reshape(subject_ids, :, 1))

write_plain(joinpath(out_dir, "beta_c.tsv"), reshape(c_beta, :, 1))
write_plain(joinpath(out_dir, "beta_d.tsv"), reshape(d_beta, :, 1))
write_plain(joinpath(out_dir, "beta_K.tsv"), K_beta)
write_plain(joinpath(out_dir, "beta_M.tsv"), M_beta)

write_plain(joinpath(out_dir, "tau_c.tsv"), reshape(c_tau, :, 1))
write_plain(joinpath(out_dir, "tau_d.tsv"), reshape(d_tau, :, 1))
write_plain(joinpath(out_dir, "tau_K.tsv"), K_tau)
write_plain(joinpath(out_dir, "tau_M.tsv"), M_tau)

write_plain(joinpath(out_dir, "beta_hat.tsv"), reshape(beta_hat, :, 1))
write_plain(joinpath(out_dir, "tau_hat.tsv"), reshape(tau_hat, :, 1))
write_plain(joinpath(out_dir, "Sigma_hat.tsv"), Sigma_hat)

# 保存维度，便于 R 端严格按维度读回矩阵
meta = DataFrame(
    m = [m],
    q_beta = [q_beta],
    q_w = [q_w],
    q_z = [q_z],
    q_sigma = [q_sigma],
    r_tau = [r_tau]
)
CSV.write(joinpath(out_dir, "meta.csv"), meta)

# 保存 η 各 nuisance direction 的标签，便于调试
eta_labels = vcat(["tau_$(j)" for j in 1:q_w], Sigma_labels)
CSV.write(joinpath(out_dir, "eta_labels.csv"), DataFrame(parameter = eta_labels))

println("Null quantities written to: ", out_dir)
