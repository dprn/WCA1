const Kernel = Array{OffsetArrays.OffsetArray{Float64,2,Array{Float64,2}},2}

# Kolmogorov kernel
K(t,p,q; ν = 1.) = sqrt(3)*exp(-4*(3*(p[1]-q[1])^2-3*(p[1]-q[1])*(p[2]+q[2])*t+t^2*(p[2]^2+p[2]*q[2]+q[2]^2))/t^3/ν^2)/π/t^2/ν^2

# Precomputation of the kernel
function kernel_computation(freq::Union{Frequencies, FloatRange{Float64}}, slopes::FloatRange{Float64}, τ;  n::Int = 20, args...)
    [centered([K(τ, (freq[k],slopes[l]), (freq[k] + i*step(freq),slopes[l] + j*step(slopes)); args...) for i in -n:n, j in -n:n]) 
            for k in 1:length(freq), l in 1:length(slopes)]
end

# Given a precomputed kernel K and a signal a, compute the required integral
function apply_kernel(a, K::Kernel)
    @assert size(K) == size(a)
    b = similar(a)
    for k in 1:size(K,1), l in 1:size(K,2)
        el = 0
        KK = K[k,l]
        for i in axes(KK,1), j in axes(KK,2)
            el += a[ clamp(k+i, 1, size(a,1)), clamp(l+j, 1, size(a,2))]*KK[i,j]
        end
        b[k,l] = el
    end
    b
end

arg(z) = atan(imag(z), real(z))
σ(z) = exp(im*arg(z))*min(1, max(abs(z), -1))

function wc_delay(input::Lift, α, β, γ; 
        K::Union{Nothing, Kernel} = nothing, τdx = 20, 
        normalize_freq::Bool = true, normalize_sigmoid::Bool = true,
        args...)
    x = time(input)
    y = normalize_freq ? normalize(freq(input)) : freq(input)
    z = slopes(input)
    
    # adapted sigmoid
    M = normalize_sigmoid ? maximum(abs.(input[:,:,:])) : 1
    σ(z) = exp(im*arg(z))*min(M, max(abs(z), -M))
    
    τ=τdx*step(x)
    if K == nothing 
        K = kernel_computation( y, z, τ; args...)
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