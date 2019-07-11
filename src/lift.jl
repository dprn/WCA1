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

function gradient_uneven(F::Array{T,2},dx::Real ,dy::Real) where T<:Number
    gradx = (F - circshift(F,[1,0]))/dx
    grady = (circshift(F,[0,-1]) - circshift(F,[0,1]))/(2*dy)
    return gradx, grady
end

# Computation of the slopes for a given input matrix M
function slopes(M::Array{T, 2}, zs::FloatRange{Float64}, 
				d_time ::Real, d_freq::Real;
                epsGrad = 1e-4, σ = (1., 1.)) where {T<:Real}
    # Recall that the xs correspond to frequencies 
    # and the ys to the time
    g_freq, g_time = gradient_uneven(M, d_freq, d_time) 
    gradNorm = g_freq.^2+g_time.^2
    
    # Generate the matrix of slopes
    N = length(zs)
    slopeMatrix = similar(M, Union{Int64,Nothing})
    for i=1:size(M,1),  j=1:size(M,2)
        if gradNorm[i,j] > epsGrad^2
            if first(zs) <= (-g_time[i,j]/g_freq[i,j]) <= last(zs)
                slopeMatrix[i,j] = round(Int,(-g_time[i,j]/g_freq[i,j])/rng(zs)*(N/2-1/2)+N/2+1/2)
            elseif (-g_time[i,j]/g_freq[i,j]) > last(zs)
                slopeMatrix[i,j] = N
            else
                slopeMatrix[i,j] = 1
            end
        else
            slopeMatrix[i,j] = nothing
        end
    end
    slopeMatrix
# 
#     imfilter(slopeMatrix,ImageFiltering.Kernel.gaussian((σ[1],σ[2]*dx/dy)))
end


# There is a problem: when computing the gradients since 
# the range of the frequencies is humungous w.r.t. the 
# time. To fix this, we compute it as if they where in the 
# range [0,1].
rng(t::Union{FloatRange{T},Frequencies, Array{T,1}}) where T<: Real = last(t) - first(t)
normalize(x :: Union{Frequencies, FloatRange{T}}) where T = (first(x)/last(x)):(step(x)/last(x)):1

function zs(m::STFT; N::Int = 100, ratio::Real = 10) 
    f_norm = normalize(freq(m))
    ratio*range(-rng(f_norm)/rng(time(m)),rng(f_norm)/rng(time(m)), length = N)
end

slopes(m::STFT; args...) = slopes(abs.(vals(m)), zs(m; args...), step(time(m)), step(normalize(freq(m))))

function lift(m::STFT; N::Int = 100, args...)
    slopeMatrix = slopes(m, N=N; args...)
    
    imgLift = zeros(eltype(m),(size(m,1),size(m,2),N))
    for i=1:size(m,1), j=1:size(m,2)
        if slopeMatrix[i,j] != nothing
            imgLift[i,j,slopeMatrix[i,j]] = m[i,j]
        else
            imgLift[i,j,:] = ones(N)*m[i,j] / N
        end
    end

    Lift(imgLift, freq(m), time(m), zs(m, N=N; args...), width(m))
end

function project(Φ::Lift)
    f = sum(Φ[:,:,:], dims = 3) |> x->dropdims(x, dims = 3)
    STFT( f, freq(Φ), time(Φ), width(Φ) )
end