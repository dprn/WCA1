export KernParams, Kern, tol, τ, b, params

# A container for the parameters of the kernel
struct KernParams
    τ :: Float64
    b :: Float64
    tol :: Float64    
end

tol(p::KernParams) = p.tol
τ(p::KernParams) = p.τ
b(p::KernParams) = p.b
params(p::KernParams) = [p.τ, p.b, p.tol]

# Explicit expression of the Kolmogorov equation kernel
function K(t,x,y; b = 1.) 
    g = 3*(x[1]-y[1])^2-3*(x[1]-y[1])*(x[2]+y[2])*t+t^2*(x[2]^2+x[2]*y[2]+y[2]^2)
    sqrt(3)/(2π*t^2*b)*exp(-g/(t^3*b))
end
K(p :: KernParams, x, y) = K(τ(p), x, y, b = b(p) )

"""
	Kern 

Container for pre-computed Kolmogorov kernel. It is a 2d `AbstractArray`,
whose value at (i,j) is a list of (I, v), where v is the value of the 
kernel at indexes I. 
"""
struct Kern <: AbstractArray{Array{Tuple{CartesianIndex{2},Float64},1}, 2}
    freq :: Union{Frequencies{Float64}, FloatRange{Float64}}
    slopes ::FloatRange{Float64}
    params :: KernParams
    vals :: Array{ Array{Tuple{CartesianIndex{2},Float64},1}, 2 }
end

Kern(f, s, τ, b, tol, vals, idx) = Kern(f,s, KernParams(τ, b, tol), vals, idx)

tol(k::Kern) = tol(k.params)
τ(k::Kern) = τ(k.params)
b(k::Kern) = b(k.params)

f(k::Kern) = k.freq
s(k::Kern) = k.slopes
params(k::Kern) = k.params
# vals(k::Kern) = k.vals

Base.size(k::Kern) = (length(k.freq), length(k.slopes))
Base.IndexStyle(::Type{<:Kern}) = IndexCartesian()
Base.getindex(k::Kern, i::Int, j::Int) = k.vals[i,j]

###############
#
#	CONSTRUCTOR
#
###############

c(τ, b, tol) = -4b*τ*(log(2π*b*τ^2/sqrt(3))+log(tol))
c(p::KernParams) = c(params(p)...)

"""
	threshold_idx(freqs, slopes, i, j, kernel_parameters)

Returns the list of positions (as an `Array{CartesianIndex{2},1}` containing
the indices w.r.t. the grid `freqs x slopes`) at which the Kolmogorov kernel
`k(freqs[i], slopes[j], *,*)` is bigger than the tolerance.

"""
#=
function threshold_idx(f, s, i, j, p) 
    ν = s[j]
	idxs = CartesianIndex{2}[]
    
    m = floor(Int, sqrt(c(p))/step(s))
    for k in max(-m, 1-j):min(m, length(s)-j)
        a = max(c(p)- (k*step(s))^2 , 0 ) # to avoid rounding errors, which make a<0
        μ_plus = τ(p)/(2*step(f))*( 2ν + k*step(s) + sqrt(a/3) )
        μ_minus = τ(p)/(2*step(f))*( 2ν + k*step(s) - sqrt(a/3) )
        C = CartesianIndex{2}[ CartesianIndex(i+ r, j+k ) for r in max( ceil(Int, μ_minus), 1-i ):min( floor(Int, μ_plus), length(f)-i)]
        append!(idxs,C)
    end
    
    idxs
end
=#

function threshold_idx(f, s, i, j, p) 
    ν = s[j]
	idxs = CartesianIndex{2}[]
    
    m = floor(Int, sqrt(c(p))/step(s))
    for k in max(-m, 1-j):min(m, length(s)-j)
        a = max(c(p)- (k*step(s))^2 , 0 ) # to avoid rounding errors, which make a<0
        μ_plus = τ(p)/(2*step(f))*(- 2ν - k*step(s) + sqrt(a/3) )
        μ_minus = τ(p)/(2*step(f))*(- 2ν - k*step(s) - sqrt(a/3) )
        C = CartesianIndex{2}[ CartesianIndex(i+ r, j+k ) for r in max( ceil(Int, μ_minus), 1-i ):min( floor(Int, μ_plus), length(f)-i)]
        append!(idxs,C)
    end
    
    idxs
end 


"""
	Kern(freqs, slopes, kernel_parameters)

Creates an instance of Kernel corresponding to the given parameters. In 
particular, it only computes it for values bigger than the tolerance.
"""
function Kern(f::Union{Frequencies, FloatRange{Float64}}, s::FloatRange{Float64}, p = KernParams(.1, 1, 1e-3))
    v =  Array{ Array{Tuple{CartesianIndex{2},Float64},1}, 2 }(undef, (length(f), length(s)))
    
    for i in 1:length(f), j in 1:length(s)
        idx = threshold_idx(f, s, i, j, p)
        v[i,j] = Tuple{CartesianIndex{2},Float64}[ (I, K(p, (f[i], s[j]), (f[I[1]], s[I[2]])) ) for I in idx ] 
    end    
    
    Kern(f,s,p,v)
end



###############
#
#	APPLICATION
#
###############

"""
	apply_kernel(a, kern)

Filters the given signal `a` via the kernel `kern`, which is of type `Kern`
"""
function apply_kernel(a, k::Kern)
	@assert size(a) == size(k)
    b = similar(a)
    for i in 1:size(a, 1), j in 1:size(a, 2)
        el = 0
        for x in k[i,j]
            el += a[ x[1] ] * x[2]
        end
        b[i,j] = el
    end
    b
end
