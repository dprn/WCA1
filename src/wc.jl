arg(z) = atan(imag(z), real(z))
σ(z) = exp(im*arg(z))*min(1, max(abs(z), -1))

function wc_delay(input::Lift, α, β, γ; 
        K::Union{Nothing, Kern} = nothing, τdx = 20, 
        normalize_freq::Bool = true, normalize_sigmoid::Bool = true,
        b = 1., tol = 1e-3,
        args...)
    x = time(input)
    y = normalize_freq ? normalize(freq(input)) : freq(input)
    z = slopes(input)
    
    # adapted sigmoid
    M = normalize_sigmoid ? maximum(abs.(input[:,:,:])) : 1
    σ(z) = exp(im*arg(z))*min(M, max(abs(z), -M))
    
    τ=τdx*step(x)
    if K == nothing 
        K = Kern( y, z, KernParams(τ, b, tol))
    elseif params(K).τ != τ
        println("WARNING: Computed τ differs from the τ of the kernel")
    end
    
    iteration(Φt, Φt_τ, h) =  (1-α*step(x))*Φt + γ*step(x)*step(y)*step(z)*apply_kernel(σ.(Φt_τ),K) + β*step(x)*h

    out = Array{Complex{Float64}, 3}(undef, length(y), length(x), length(z))
    out[:,1,:] = copy(input[:, 1,:])
    mat_zero = zeros(Complex{Float64},size(out[:,1,:])...)

    p = Progress(length(x), desc = "WC evolution...")

    for t in 1:τdx
         out[:, t+1,:] = iteration(out[:,t,:],mat_zero, input[:,t,:])

         next!(p)
    end
    for t in τdx+1:(length(x)-1)
         out[:,t+1,:] = iteration(out[:,t,:],out[:,t-τdx,:], input[:,t,:])

         next!(p)
    end
    Lift(out, freq(input), x, z, width(input))
end