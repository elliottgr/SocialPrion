## Model parameters
	
	r = 1.0
	G = 0.5
	g = 0.5
	K_B = 1.0
	K_y = 1.0
	K_G = 1.0
	H = 0.1
	B = 1.0
	p = 0.2
	E = 0.1
	
	## Plotting parameters

	N_min = 0.0
	N_step = 0.01
	N_max = 1.0


	## Base model
	
	function dBdt(r, B, G, g, K_B, K_y, K_G, H, p, E)
		return r*B*(1-(B/K_B)) - (B*G*E)
	end
	
	function dGdt(r, B, G, g, K_B, K_y, K_G, H, p, E)
		return r*G*(1-((G+g)/K_y)) + (H * g) - (B*p*G) 
	end
	
	function dgdt(r, B, G, g, K_B, K_y, K_G, H, p, E)
		return r*g*(1-((G+g)/K_y)) + (B*p*G) - (H*g)
	end

	## Glucose limited model

	function dBdt_G(r, B, G, g, K_B, K_y, K_G, H, p, E)
		return r*B*(1-(B/K_B)) - (B*G*E)
	end
	
	function dGdt_G(r, B, G, g, K_B, K_y, K_G, H, p, E)
		return r*G*(1 - ((G+B)/K_G) - ((g+G)/K_y)) + (H * g) - (B*p*G) 
	end

	function dgdt_G(r, B, G, g, K_B, K_y, K_G, H, p, E)
		return r*g*(1 - ((g+G)/K_y)) + (B*p*G) - (H*g)
	end
	
	## Plotting
	p1 = plot(title = "Base Model", ylim = (-0.5, 0.5), legend = :none)
	p2 = plot(title = "Glucose Limited Model", ylim = (-0.5, 0.5), yformatter = :none)
	
	plot!(p1, N_min:N_step:N_max, [dBdt(r, N, G, g, K_B, K_y, K_G, H, p, E) for N in N_min:N_step:N_max], label = "B")
	
	plot!(p1, N_min:N_step:N_max, [dGdt(r, B, N, g, K_B, K_y, K_G, H, p, E) for N in N_min:N_step:N_max], label = "G")

	plot!(p1, N_min:N_step:N_max, [dgdt(r, B, G, N, K_B, K_y, K_G, H, p, E) for N in N_min:N_step:N_max], label = "g")

	plot!(p2, N_min:N_step:N_max, [dBdt_G(r, N, G, g, K_B, K_y, K_G, H, p, E) for N in N_min:N_step:N_max], label = "B")

	plot!(p2, N_min:N_step:N_max, [dGdt_G(r, B, N, g, K_B, K_y, K_G, H, p, E) for N in N_min:N_step:N_max], label = "G")

	plot!(p2, N_min:N_step:N_max, [dgdt_G(r, B, G, N, K_B, K_y, K_G, H, p, E) for N in N_min:N_step:N_max], label = "g")
	
	plot(p1, p2, layout = (1,2))


	# plot([base_B, base_G, base_g, glucose_B, glucose_G, glucose_G], layout = (2,1))
