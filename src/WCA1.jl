module WCA1

using WAV, FFTW, OffsetArrays, Plots

import DSP, Statistics
using DSP.Windows, ImageFiltering
import FFTW: Frequencies, fftfreq, rfftfreq

export stft, istft, show_stft, lift, project, wc_delay,
	Lift, STFT, kernel_computation, freq, slopes, time, normalize

include("stft.jl")
include("lift.jl")
include("wc.jl")

end # module
