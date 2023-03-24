using DifferentialEquations, Plots

include("prionStructs.jl")


## version of function that uses vectors
# function f(du, u, p ,t)
#     B_t, G_t, g_t = u
#     b_B, b_g, b_G, d_B, d_g, d_G, E, H, pr = p ##i hate this! 
#     du[1] = b_B*B_t - G_t*E - B_t*d_B
#     du[2] = b_G*G_t + H*g_t - B_t*pr*G_t - d_G*G_t
#     du[3] = b_g*g_t + B_t*pr*G_t - H*g_t - d_g*g_t
# end

## rewritten to use structs ()
# function f(du, u, p ,t)
    B_t, G_t, g_t = u.pop
    # b_B, b_g, b_G, d_B, d_g, d_G, E, H, pr = p
    du[1] = p.b_B*u.B_t - u.G_t*p.E - u.B_t*p.d_B
    du[2] = p.b_G*u.G_t + p.H*u.g_t - u.B_t*p.pr*u.G_t - p.d_G*u.G_t
    du[3] = p.b_g*u.g_t + u.B_t*p.pr*u.G_t - p.H*u.g_t - p.d_g*u.g_t
# end


# function f(du, u, p ,t)
#     # B_t, G_t, g_t = u
#     # b_B, b_g, b_G, d_B, d_g, d_G, E, H, pr = p
#     du[1] = p.b_B*u.B_t - u.G_t*p.E - u.B_t*p.d_B
#     du[2] = p.b_G*u.G_t + p.H*u.g_t - u.B_t*p.pr*u.G_t - p.d_G*u.G_t
#     du[3] = p.b_g*u.g_t + u.B_t*p.pr*u.G_t - p.H*u.g_t - p.d_g*u.g_t
# end

function main()
    ## Parameters
    b_B = 1.0 ## bacterial growth rate
    b_G = 1.0 ## GRA+ growth rate
    b_g = 1.0 ## gra- growth rate

    d_B = 1.0
    d_G = 1.0
    d_g = 1.0

    E = 0.1 ## Effect of ethanol on bacterial death
    H = 0.01 ## effect of HSP on prion + yeast
    pr = 0.05 ## effect of bacteria on prion activation

    p = [b_B, b_g, b_G, d_B, d_g, d_G, E, H, pr]

    B_0 = 100 ## Bacterial pop size // Guessing at a reasonable population size
    G_0 = 100 ## GRA+ pop size
    g_0 = 1 ## gra- pop size

    u = [B_0, G_0, g_0]
    # u = state_var(B_0, G_0, g_0)

    tspan = (0.0, 3.0)
    prob = ODEProblem(f, u, tspan, p)
    sol = solve(prob)

    plot(sol, labels = ["B(t)" "G(t)" "g(t)"])


end