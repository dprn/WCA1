struct Lift{T}
    lift::Array{T, 3}
    freq::Frequencies
    time::FloatRange{Float64}
    slopes::FloatRange{Float64}
    width::Int
end

Base.eltype(m::Lift{T}) where T = eltype(lift(m))
Base.size(m::Lift{T}, etc...) where T = size(lift(m), etc...)
Base.getindex(m::Lift{T}, etc...) where T = getindex(lift(m), etc...)

width(L::Lift) = L.width
freq(L::Lift) = L.freq
lift(L::Lift) = L.lift
Libc.time(L::Lift) = L.time
slopes(L::Lift) = L.slopes

# function grad(I)
#     g1, g2 = KernelFactors.ando3()
#     imfilter(I, g1), imfilter(I, g2)
# end

grad(I) = imgradients(I, KernelFactors.ando3)

function compute_slopes(SS; threshold = 1e-3, args...)
    M = abs.(vals(SS))
    gx, gy = grad(M)
    
    gx = gx * length(SS.freq)
    gy = gy / step(SS.time)
    
    G = similar(M)
    for i in 1:length(M)
        G[i] = abs(gx[i]) > threshold  ? -gy[i]/gx[i] : 0.
    end
    G
end

# Functions to cut away too small gradients 
function cut(G, σ1, σ2 = σ1) 
    M = similar(G)
    
    m = Statistics.mean(G)
    s = Statistics.std(G)
    @simd for i in 1:length(M)
        if m-σ1*s < G[i] < m+σ2*s
            M[i] = G[i]
        else
            M[i] = 0.
        end
    end
    M
end

function auto_cut(G; p = 0.95, args...)
    n(M) = sum(abs.(M-G))/sum(abs.(G))
    C = ([ n(cut(G, σ1, σ2)) for σ1 in 1:10, σ2 in 1:10 ]).-(1-p)
    cut(G,Tuple(findmin(abs.(C))[2])...)
end


normalize(x) = (first(x)/last(x)):(step(x)/last(x)):1

rng(t) = last(t) - first(t)

zs(M, νMin, νMax, N) = range(νMin,νMax, length = N)

function compute_slope_matrix(M, νMin, νMax, N = 100; args...)
    Z = zs(M, νMin, νMax, N)
    slopeMatrix = similar(M, Union{Int,Nothing})
    A = round.(Int,(M .- first(Z))*(N-1)/(rng(Z)) .+ 1)
    for i in 1:length(M)
        if first(Z) <= M[i] <= last(Z)
            slopeMatrix[i] = A[i]
        elseif M[i] > last(Z)
                slopeMatrix[i] = nothing#N
        elseif M[i] < first(Z)
                    slopeMatrix[i] = nothing#1
        else
            slopeMatrix[i] = nothing
        end
    end
    slopeMatrix, Z
end


function slopes(SS::STFT, νMin, νMax, N = 100; args...) where {T<:Real}
    G = compute_slopes(SS; args...)
    #M = auto_cut(G; args...)                    ### Placeholder for automated slope detection
    compute_slope_matrix(G, νMin, νMax, N)
end

function lift(m::STFT; νMin=-1, νMax=1, N::Int = 100,  args...)
    slopeMatrix, Z = slopes(m, νMin, νMax, N; args...)

    imgLift = zeros(eltype(m),(size(m,1),size(m,2),N))
    for i=1:size(m,1), j=1:size(m,2)
        if slopeMatrix[i,j] != nothing
            imgLift[i,j,slopeMatrix[i,j]] = m[i,j]
        else
            imgLift[i,j,:] = ones(N)*m[i,j] / N
        end
    end

    Lift(imgLift, freq(m), time(m), Z, width(m))
end

function project(Φ::Lift)
    f = sum(Φ[:,:,:], dims = 3) |> x->dropdims(x, dims = 3)
    STFT( f, freq(Φ), time(Φ), width(Φ) )
end

show_stft(m::Lift; args...) = show_stft(project(m); args...)
show_istft(m::Lift; args...) = show_istft(project(m); args...)
