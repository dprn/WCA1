
# Input sound

rate = 16000.
duration = 2
samples = round(Int,duration*rate)

x = [ samples/7 <= t <= 8*samples/14 ? sin(3000*2*π*t/rate+2pi*150*sin(2*pi*t/rate)) : 0. for t in 1:samples ]

## Short time Fourier transform

m =stft(x, 500, 450, fs = rate, window = hanning)


## smoothing

using ImageFiltering

λ=2

M = STFT(0.01.*ComplexF64.(imfilter(abs.(m.stft),ImageFiltering.Kernel.gaussian((λ,2λ)))), m.freq, m.time, m.width)
                                                    

## Lift

Lm = lift(M; νMin=-1.5, νMax=1.5, N=100)

## WC evolution

χ = 20

α = 53
β = 1
γ = 55

τ = χ * step(time(Lm))

b = .2

k = Kern(normalize(freq(Lm)), slopes(Lm), KernParams(τ, b, 1e-6));

W = wc_delay(Lm, α, β, γ, K=k,τdx = χ) |> project

## Save results

try
    mkpath("nonlinear-chirp-results")
catch
end

cd("nonlinear-chirp-results")

save_result(M, W, α, β, γ, χ)

cd("..")

