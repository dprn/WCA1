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

function compute_slopes(SS; threshold = 1e-3)
    M = abs.(vals(SS))
    gx, gy = grad(M)
    
#     gx = gx / step(SS.freq)
    gx = gx * length(SS.freq)
    gy = gy / step(SS.time)
    
    G = similar(M)
    for i in 1:length(M)
        G[i] = gx[i] > threshold  ? -gy[i]/gx[i] : 0.
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

zs(M, N; k = 3) = range(k*minimum(M),k*maximum(M), length = N)

function compute_slope_matrix(M, N = 100)
    Z = zs(M, N)
    slopeMatrix = similar(M, Union{Int,Nothing})
    A = round.(Int,(M .- first(Z))*(N-1)/(rng(Z)) .+ 1)
    for i in 1:length(M)
        if first(Z) <= M[i] <= last(Z)
            slopeMatrix[i] = A[i]
        elseif M[i] > last(Z)
                slopeMatrix[i] = N
        elseif M[i] < first(Z)
                    slopeMatrix[i] = 1
        else
            slopeMatrix[i] = nothing
        end
    end
    slopeMatrix, Z
end

function slopes(SS::STFT, N = 100; args...) where {T<:Real}
    G = compute_slopes(SS; args...)
    M = auto_cut(G; args...)
    compute_slope_matrix(M, N)
end

function lift(m::STFT; N::Int = 100, args...)
    slopeMatrix, Z = slopes(m, N; args...)

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

# function gradient_uneven(F::Array{T,2},dx::Real ,dy::Real) where T<:Number
#     gradx = (F - circshift(F,[1,0]))/dx
#     grady = (circshift(F,[0,-1]) - circshift(F,[0,1]))/(2*dy)
#     return gradx, grady
# end

# # Computation of the slopes for a given input matrix M
# function slopes(M::Array{T, 2}, zs::FloatRange{Float64}, 
# 				d_time ::Real, d_freq::Real;
#                 epsGrad = 1e-4, σ = (1., 1.)) where {T<:Real}
#     # Recall that the xs correspond to frequencies 
#     # and the ys to the time
#     g_freq, g_time = gradient_uneven(M, d_freq, d_time) 
#     gradNorm = g_freq.^2+g_time.^2
    
#     # Generate the matrix of slopes
#     N = length(zs)
#     slopeMatrix = similar(M, Union{Int64,Nothing})
#     for i=1:size(M,1),  j=1:size(M,2)
#         if gradNorm[i,j] > epsGrad^2
#             if first(zs) <= (-g_time[i,j]/g_freq[i,j]) <= last(zs)
#                 slopeMatrix[i,j] = round(Int,(-g_time[i,j]/g_freq[i,j])/rng(zs)*(N/2-1/2)+N/2+1/2)
#             elseif (-g_time[i,j]/g_freq[i,j]) > last(zs)
#                 slopeMatrix[i,j] = N
#             else
#                 slopeMatrix[i,j] = 1
#             end
#         else
#             slopeMatrix[i,j] = nothing
#         end
#     end
#     slopeMatrix
# # 
# #     imfilter(slopeMatrix,ImageFiltering.Kernel.gaussian((σ[1],σ[2]*dx/dy)))
# end


# # There is a problem: when computing the gradients since 
# # the range of the frequencies is humungous w.r.t. the 
# # time. To fix this, we compute it as if they where in the 
# # range [0,1].
# rng(t::Union{FloatRange{T},Frequencies, Array{T,1}}) where T<: Real = last(t) - first(t)
# normalize(x :: Union{Frequencies, FloatRange{T}}) where T = (first(x)/last(x)):(step(x)/last(x)):1

# function zs(m::STFT; N::Int = 100, ratio::Real = 10) 
#     f_norm = normalize(freq(m))
#     ratio*range(-rng(f_norm)/rng(time(m)),rng(f_norm)/rng(time(m)), length = N)
# end

# slopes(m::STFT; args...) = slopes(abs.(vals(m)), zs(m; args...), step(time(m)), step(normalize(freq(m))))

# function lift(m::STFT; N::Int = 100, args...)
#     slopeMatrix = slopes(m, N=N; args...)
    
#     imgLift = zeros(eltype(m),(size(m,1),size(m,2),N))
#     for i=1:size(m,1), j=1:size(m,2)
#         if slopeMatrix[i,j] != nothing
#             imgLift[i,j,slopeMatrix[i,j]] = m[i,j]
#         else
#             imgLift[i,j,:] = ones(N)*m[i,j] / N
#         end
#     end

#     Lift(imgLift, freq(m), time(m), zs(m, N=N; args...), width(m))
# end

function project(Φ::Lift)
    f = sum(Φ[:,:,:], dims = 3) |> x->dropdims(x, dims = 3)
    STFT( f, freq(Φ), time(Φ), width(Φ) )
end