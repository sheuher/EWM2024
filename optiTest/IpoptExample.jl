using Ipopt
using Test

# hs071
# min x1 * x4 * (x1 + x2 + x3) + x3
# st  x1 * x2 * x3 * x4 >= 25
#     x1^2 + x2^2 + x3^2 + x4^2 = 40
#     1 <= x1, x2, x3, x4 <= 5
# Start at (1,5,5,1)
# End at (1.000..., 4.743..., 3.821..., 1.379...)

function eval_f(x::Vector{Float64})
    return x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
end

function eval_g(x::Vector{Float64}, g::Vector{Float64})
    g[1] = x[1] * x[2] * x[3] * x[4]
    return g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
end

function eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64})
    grad_f[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3])
    grad_f[2] = x[1] * x[4]
    grad_f[3] = x[1] * x[4] + 1
    return grad_f[4] = x[1] * (x[1] + x[2] + x[3])
end

function eval_jac_g(
    x::Vector{Float64}, 
    rows::Vector{Int32},
    cols::Vector{Int32},
    values::Union{Nothing,Vector{Float64}}
)
    if values === nothing
        # Constraint (row) 1
        rows[1] = 1
        cols[1] = 1
        rows[2] = 1
        cols[2] = 2
        rows[3] = 1
        cols[3] = 3
        rows[4] = 1
        cols[4] = 4
        # Constraint (row) 2
        rows[5] = 2
        cols[5] = 1
        rows[6] = 2
        cols[6] = 2
        rows[7] = 2
        cols[7] = 3
        rows[8] = 2
        cols[8] = 4
    else
         # Constraint (row) 1
         values[1] = x[2] * x[3] * x[4]  # 1,1
         values[2] = x[1] * x[3] * x[4]  # 1,2
         values[3] = x[1] * x[2] * x[4]  # 1,3
         values[4] = x[1] * x[2] * x[3]  # 1,4
         # Constraint (row) 2
         values[5] = 2 * x[1]  # 2,1
         values[6] = 2 * x[2]  # 2,2
         values[7] = 2 * x[3]  # 2,3
         values[8] = 2 * x[4]  # 2,4
    end
    return
end

function eval_h(
    x::Vector{Float64},
    rows::Vector{Int32},
    cols::Vector{Int32},
    obj_factor::Float64,
    lambda::Vector{Float64},
    values::Union{Nothing,Vector{Float64}}
)
    if values === nothing
        # Symmetric matrix, fill the lower left triangle only
        idx = 1
        for row in 1:4
            for col in 1:row
                rows[idx] = row
                cols[idx] = col
                idx += 1
            end
        end
    else
        # Again, only lower left triangle
        # Objective
        values[1] = obj_factor * (2 * x[4])  # 1,1
        values[2] = obj_factor * (x[4])  # 2,1
        values[3] = 0                      # 2,2
        values[4] = obj_factor * (x[4])  # 3,1
        values[5] = 0                      # 3,2
        values[6] = 0                      # 3,3
        values[7] = obj_factor * (2 * x[1] + x[2] + x[3])  # 4,1
        values[8] = obj_factor * (x[1])  # 4,2
        values[9] = obj_factor * (x[1])  # 4,3
        values[10] = 0                     # 4,4
        # First constraint
        values[2] += lambda[1] * (x[3] * x[4])  # 2,1
        values[4] += lambda[1] * (x[2] * x[4])  # 3,1
        values[5] += lambda[1] * (x[1] * x[4])  # 3,2
        values[7] += lambda[1] * (x[2] * x[3])  # 4,1
        values[8] += lambda[1] * (x[1] * x[3])  # 4,2
        values[9] += lambda[1] * (x[1] * x[2])  # 4,3

        # Second constraint
        values[1] += lambda[2] * 2  # 1,1
        values[3] += lambda[2] * 2  # 2,2
        values[6] += lambda[2] * 2  # 3,3
        values[10] += lambda[2] * 2  # 4,4
    end
    return 
end

n = 4
x_L = [1., 1., 1., 1.]
x_U = [5., 5., 5., 5.]

m = 2
g_L = [25., 40.]
g_U = [2.0e19, 40.]

prob = Ipopt.CreateIpoptProblem(
    n,
    x_L,
    x_U,
    m,
    g_L,
    g_U,
    8,
    10,
    eval_f,
    eval_g,
    eval_grad_f,
    eval_jac_g,
    eval_h,
)

prob.x = [1.0, 5.0, 5.0, 1.0]

#tatus = Ipopt.IpoptSolve(prob)

function intermediate(
    alg_mod::Cint,
    iter_count::Cint,
    obj_value::Float64,
    inf_pr::Float64,
    inf_du::Float64,
    mu::Float64,
    d_norm::Float64,
    regularization_size::Float64,
    alpha_du::Float64,
    alpha_pr::Float64,
    ls_trials::Cint,
)
    m, n = 2, 4
    x, z_L, z_U = zeros(n), zeros(n), zeros(n)
    g, lambda = zeros(m), zeros(m)
    scaled = false
    Ipopt.GetIpoptCurrentIterate(prob, scaled, n, x, z_L, z_U, m, g, lambda)
    x_L_violation, x_U_violation = zeros(n), zeros(n)
    compl_x_L, compl_x_U, grad_lag_x = zeros(n), zeros(n), zeros(n)
    nlp_constraint_violation, compl_g = zeros(m), zeros(m)
    Ipopt.GetIpoptCurrentViolations(
        prob,
        scaled,
        n,
        x_L_violation,
        x_U_violation,
        compl_x_L,
        compl_x_U,
        grad_lag_x,
        m,
        nlp_constraint_violation,
        compl_g,
    )
    @test x .+ x_L_violation >= x_L
    @test x .- x_U_violation <= x_U
    @test g_L <= g .- nlp_constraint_violation <= g_U
    return iter_count < 1  # Interrupts after one iteration.
end

#Ipopt.SetIntermediateCallback(prob, intermediate)

solvestat = Ipopt.IpoptSolve(prob)


#println(Ipopt.ApplicationReturnStatus[status])
println(prob.x)
println(prob.obj_val)

prob.x
a=zeros
eval_g(prob.x, zeros(4))