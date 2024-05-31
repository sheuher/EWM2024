using Ipopt


# min x1^2 - 1/2*x1 - x2 - 2
# st  x1^2 - 4x1 + x2 + 1 <= 0
#     1/2x1^2 + x2^2 - x1 - 4 <= 0


function eval_f(x)
    return x[1]^2 - .5*x[1] - x[2] - 2
end

function eval_g(x, g)
    g[1] = x[1]^2 - 4*x[1] + x[2] + 1 
    return g[2] = .5*x[1]^2 + x[2]^2 - x[1] - 4
end

function eval_grad_f(x, grad_f)
    grad_f[1] = 2*x[1] - 0.5
    return grad_f[2] = -1
end

function eval_jac_g(
    x, 
    rows,
    cols,
    values
)
    if values === nothing
        # m = 2
        # n = 2
        rows[1] = 1
        cols[1] = 1
        rows[2] = 2
        cols[2] = 1
        rows[3] = 1
        cols[3] = 2
        rows[4] = 2
        cols[4] = 2
    else
        values[1] = 2*x[1] - 4
        values[2] = x[1] - 1
        values[3] = 1
        values[4] = 2*x[2]
    end
    return
end

nzJ = 4
nzH = 0
n = 2
m = 2
x_L = -5*ones(n)
x_U = 5*ones(n)
g_L = -Inf*ones(n)
g_U = zeros(n)

prob = Ipopt.CreateIpoptProblem(
    n,
    x_L,
    x_U,
    m,
    g_L,
    g_U,
    nzJ,
    nzH,
    eval_f,
    eval_g,
    eval_grad_f,
    eval_jac_g,
    nothing
)

AddIpoptStrOption(prob, "hessian_approximation", "limited-memory")

prob.x = ones(n)

IpoptSolve(prob)

prob.x
"""
2-element Vector{Float64}:
 1.0623762746384762
 2.120861759683611"""
prob.obj_val