using DifferentialEquations

## Parameters
b_B = 1.1 ## bacterial growth rate
b_g = 1.1 ## gra- growth rate
b_G = 1.1 ## GRA+ growth rate
d_B = .9
d_g = .9
d_G = .9
E = 0.05 ## Effect of ethanol on bacterial death
H = 0.01 ## effect of HSP on prion + yeast
p = 0.05 ## effect of bacteria on prion activation

B_0 = 10^9 ## Bacterial pop size // Guessing at a reasonable population size
g_0 = 0 ## gra- pop size
G_0 = 10^6 ## GRA+ pop size


ode = @ode_def GRAModel begin
    dg = b_g*g_t + B_t*p*G_t - H*g_t - d_g*g_t
    dG = b_G*G_t + H*g_t - B_t*p*G_t - d_G*G_t
    dB = b_B*B_t - B_t*E*d_B
end
