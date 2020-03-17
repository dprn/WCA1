using Dates, Images

function save_result(m, α, β, γ, τdx; rate = 16000.) 
	str = string(Dates.format(now(),"yyyy:mm:dd"), " at ", Dates.format(now(),"HH"),"h",Dates.format(now(),"MM")," - α=",α,", β=",β,", γ=",γ,", τdx=",τdx)

	save(string(str," - STFT.png"), show_stft(m))
	save(string(str," - Wave.png"), show_istft(m))

	sound = istft(project(m))
	sound = 2*(sound .- minimum(sound))/(maximum(sound)-minimum(sound)).-1
	wavwrite(sound, string(str," - sound.wav"), Fs = rate)
end

using WAV

import Pkg
Pkg.activate("../../WCA1")
using WCA1

# include("linear_chirp.jl")
include("bars.jl")