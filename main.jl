using DifferentialEquations, Plots

## version of function that uses vectors
# function f(du, u, p ,t)
#     B_t, G_t, g_t = u
#     b_B, b_g, b_G, d_B, d_g, d_G, E, H, pr = p ##i hate this! 
#     du[1] = b_B*B_t - G_t*E - B_t*d_B
#     du[2] = b_G*G_t + H*g_t - B_t*pr*G_t - d_G*G_t
#     du[3] = b_g*g_t + B_t*pr*G_t - H*g_t - d_g*g_t
# end

## version with population maxes
function f(du, u, p ,t)
    B_t, G_t, g_t = u
    r_B, r_g, r_G, k_B, k_Y, E, H, pr = p ##i hate this! 
    du[1] = r_B*B_t*(1-B_t/k_B) - B_t*G_t*E # dB/dT
    du[2] = r_G*G_t*(1-((G_t+g_t)/k_Y)) + H*g_t*B_t - B_t*pr*G_t # dG/dt
    du[3] = r_g*g_t*(1-((G_t+g_t)/k_Y)) + B_t*pr*G_t - H*g_t*B_t # dg/dT
end



function main()
    ## Parameters
    
    r_B = 1.0
    r_g = 1.1
    r_G = 1.1

    k_B = 1.0
    k_Y = 1.0

    E = 0.1 ## Effect of ethanol on bacterial death
    H = 0.01 ## effect of HSP on prion + yeast
    pr = 0.05 ## effect of bacteria on prion activation

    p = [r_B, r_g, r_G, k_B, k_Y, E, H, pr]

    B_0 = 1 ## Bacterial pop size // Guessing at a reasonable population size
    G_0 = 0.1 ## GRA+ pop size
    g_0 = 0.01 ## gra- pop size

    u = [B_0, G_0, g_0]
    # u = state_var(B_0, G_0, g_0)

    tspan = (0.0, 3.0)
    prob = ODEProblem(f, u, tspan, p)
    sol = solve(prob)

    plot(sol, labels = ["B(t)" "G(t)" "g(t)"])


end