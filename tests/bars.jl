# Simple bars for tests

using ImageFiltering

rate = 16000.
duration = 2
samples = round(Int,duration*rate)
y = [ t <= 8*samples/14 ? sin(1000*2*π*t/rate+1000*t/rate*2*π*t/rate) : 0. for t in 1:samples ]

SS = stft(y, 500, 450, fs = rate);

function two_bars(p1, p2, km1, km2, N1, N2; λ = 4)
    # the bars
    U0v1=10*[(((i+10)p1+km1-1<j<(i-10)p1+km1+1) && (N2/6<i<N2/2) ) ? 1. : 0. for j in 1:N1, i in 1:N2]
    U0v2=10*[(((i+3)p2+km2-1<j<(i-3)p2+km2+1) && (N2/6<i<N2/3) ) ? 1. : 0. for j in 1:N1, i in 1:N2]

    # a little bit of smoothing
    imfilter!(U0v1, U0v1,ImageFiltering.Kernel.gaussian((λ,3λ)))
    imfilter!(U0v2, U0v2,ImageFiltering.Kernel.gaussian((λ,3λ)))

    # combining them
    U0v1/maximum(U0v1)+U0v2/maximum(U0v2)
end

# wrapping everything up
function two_bars_stft(p1, p2, km1, km2, SS::STFT; args...) 
    U0 = two_bars(p1,p2,km1,km2,  length(SS.freq), length(SS.time); args...)
    STFT(ComplexF64.(U0), SS.freq, SS.time, SS.width)
end

p0 = 1000 .*step(SS.time)/step(SS.freq)
SS = two_bars_stft(p0/2, -2p0, 150, 100, SS)