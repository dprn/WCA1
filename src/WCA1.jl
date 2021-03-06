module WCA1

using WAV, FFTW, OffsetArrays, Plots, ProgressMeter

import DSP, Statistics
using DSP.Windows, ImageFiltering
import FFTW: Frequencies, fftfreq, rfftfreq

export stft, istft, show_stft, lift, project, wc_delay,
	Lift, STFT, freq, slopes, time, normalize, show_istft

include("stft.jl")
include("lift.jl")
include("kernel.jl")
include("wc.jl")

end # module
