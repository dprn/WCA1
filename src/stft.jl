# We will use the already defined functions in DSP, 
# adapting the spectrogram and stft functions.
# 
# First of all we define a container type for our stft.

const FloatRange{T} = StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}

struct STFT{T, F<:Union{Frequencies, AbstractRange}} 
    stft::Matrix{T}
    freq::F
    time::FloatRange{Float64}
    width::Int
end

Libc.time(s::STFT) = s.time
freq(s::STFT) = s.freq
vals(s::STFT) = s.stft
# window_norm(s::STFT) = s.window_norm
width(s::STFT) = s.width

Base.eltype(m::STFT) = eltype(vals(m))
Base.size(m::STFT{T, Frequencies}, etc...) where T = size(vals(m), etc...)
Base.getindex(m::STFT{T, Frequencies}, etc...) where T = getindex(vals(m), etc...)
Base.extrema(m::STFT) = extrema(vals(m))

# Functions to recover informations about the stft from the 
# quantities we saved inside the type
fs(s::STFT) = round(Int, width(s)/(2*first(time(s))))
noverlap(s::STFT) = round(Int, width(s) - fs(s)*step(time(s)))
or_length(m::STFT) = (length(time(m))-1)*50 + width(m)

# We can now define the stft function, which is just a 
# wrapper around DSP.stft.

function stft(s::AbstractVector{T}, n::Int=length(s)>>3, 
					noverlap::Int=n>>1;
                    onesided::Bool=eltype(s)<:Real, nfft::Int=DSP.nextfastfft(n), 
                    fs::Real=1, window::Union{Function,AbstractVector,Nothing}=nothing) where T

    out = DSP.stft(s, n, noverlap; onesided=onesided, fs=fs, window=window)
    STFT(out, onesided ? rfftfreq(nfft, fs) : fftfreq(nfft, fs),
                (n/2 : n-noverlap : (size(out,2)-1)*(n-noverlap)+n/2) / fs, n)
end

function show_stft(m::STFT, pre = x->x) 
	# Function to plot an STFT
    T = map(x->round(x, digits = 2), range(first(time(m)), last(time(m)), length= 4))
    F = map(x->round(Int, x), range(first(freq(m)), last(freq(m)), length= 5))
    heatmap(pre.(abs.(vals(m))), xaxis = ("Time (s)"), yaxis = ("Frequencies (Hz)"),
        xticks = (range(1,size(m,2), length= 4), T),
        yticks = (range(1,size(m,1), length= 5), F))
end

# We now define the inverse STFT.

function addto!(dest, dof, src, sof)
	# Function to add vectors (inspired by Compat.copyto!)
    @assert sof <= length(src) 
    @assert dof + length(src) - sof  <= length(dest)
    @inbounds for i in 0:(length(src) - sof)
        dest[i + dof] += src[i + sof]
    end
end

function istft(m::STFT)
    onesided = all(freq(m) .>= 0)
    @assert onesided
    a = irfft(vals(m), width(m), 1)
    x = zeros(or_length(m))
    tmp = zeros(or_length(m))
    for i in 1:size(a,2)
        addto!(x, (i-1)*(width(m)-noverlap(m))+1, a[:,i], 1)
        addto!(tmp, (i-1)*(width(m)-noverlap(m))+1, ones(size(a,1)), 1)
    end
    x./=(tmp)/2
end