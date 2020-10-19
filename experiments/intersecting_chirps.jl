# Input sound

rate = 16000.
duration = 2
samples = round(Int,duration*rate)

x1 = [ samples/7 <= t <= 6*samples/14 ?     sin(2000*2*π*t/rate + 1000*t/rate*2*π*t/rate)  : 0. for t in 1:samples ]
x2 = [ samples/7 <= t <= 6*samples/14 ?     sin(6000*2*π*t/rate - 1000*t/rate*2*π*t/rate)  : 0. for t in 1:samples ]


x = x1+x2

## Short time Fourier transform

m =stft(x, 500, 450, fs = rate, window = hanning)


## smoothing

using ImageFiltering

λ=2

M = STFT(0.01.*ComplexF64.(imfilter(abs.(m.stft),ImageFiltering.Kernel.gaussian((λ,2λ)))), m.freq, m.time, m.width)
                                                    

## Lift

Lm = lift(M; νMin=-1, νMax=1, N=100)

## WC evolution

χ = 20

α = 53
β = 1
γ = 55

τ = χ * step(time(Lm))

b = .01

k = Kern(normalize(freq(Lm)), slopes(Lm), KernParams(τ, b, 1e-6));

W = wc_delay(Lm, α, β, γ, K=k,τdx = χ) |> project

## Save results

try
    mkpath("intersecting-chirps-results")
catch
end

cd("intersecting-chirps-results")

save_result(M, W, α, β, γ, χ)

cd("..")
