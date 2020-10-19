
# Input sound

rate = 16000.
duration = 2
samples = round(Int,duration*rate)

x1 = [ samples/7 <= t <= 2.5*samples/7 ? sin(2000*2*π*t/rate+1000*t/rate*2*π*t/rate) : 0. for t in 1:samples ]
x2= [ 3*samples/7 <= t <= 4.5*samples/7 ? sin(2000*2*π*t/rate+1000*t/rate*2*π*t/rate) : 0. for t in 1:samples ]

x = x1+x2

## Short time Fourier transform

m =stft(x, 500, 450, fs = rate, window = hanning)


## smoothing

using ImageFiltering

λ=2

M = STFT(0.01.*ComplexF64.(imfilter(abs.(m.stft),ImageFiltering.Kernel.gaussian((λ,2λ)))), m.freq, m.time, m.width)
                                                    

## Lift

Lm = lift(M; νMin=-0.5, νMax=1.5, N=100)

## WC evolution

χ = 20

α = 55
β = 1
γ = 55

τ = χ * step(time(Lm))

b = .05

k = Kern(normalize(freq(Lm)), slopes(Lm), KernParams(τ, b, 1e-6));

W = wc_delay(Lm, α, β, γ, K=k,τdx = χ) |> project

## Save results

try
    mkpath("interrupted-chirp-results")
catch
end

cd("interrupted-chirp-results")

save_result(M, W, α, β, γ, χ)

cd("..")
