module WCA1

using WAV, FFTW, OffsetArrays, Plots

import DSP, DSP.Frequencies, ImageFiltering.centered
using DSP.Windows

export stft, istft, show_stft, lift, project, wc_delay,
		Lift, STFT, kernel_computation, freq, slopes, time, normalize

include("stft.jl")
include("lift.jl")
include("wc.jl")

end # module
