using Dates

function save_result(m_in, m_out, α, β, γ, τdx; rate = 16000.) 
	str = string(Dates.format(now(),"yyyy-mm-dd"), " at ", Dates.format(now(),"HH"),"h",Dates.format(now(),"MM")," - a=",α,", b=",β,", g=",γ,", tdx=",τdx)

	save(string(str," - input - STFT.png"), show_stft(m_in))
	save(string(str," - input - Wave.png"), show_istft(m_in))

	save(string(str," - output - STFT.png"), show_stft(m_out))
	save(string(str," - output - Wave.png"), show_istft(m_out))
    
    normalize(sound) = 2*(sound .- minimum(sound))/(maximum(sound)-minimum(sound)).-1

	sound_in = istft(m_in) |> normalize
	wavwrite(sound_in, string(str," - input - sound.wav"), Fs = rate)

	sound_out = istft(m_out) |> normalize
	wavwrite(sound_out, string(str," - output - sound.wav"), Fs = rate)
end

using WAV
using DSP.Windows

import Pkg
Pkg.activate("../../WCA1")
using WCA1

using Images

# TESTS

try
    mkpath("results")
catch
end

cd("results")

println("=== LINEAR CHIRP ===")
include("linear_chirp.jl")

println("=== INTERRUPTED CHIRP ===")
include("interrupted_chirp.jl")

println("=== INTERSECTING CHIRPS ===")
include("intersecting_chirps.jl")

println("=== NON-LINEAR CHIRP ===")
include("nonlinear_chirp.jl")

cd("..")
