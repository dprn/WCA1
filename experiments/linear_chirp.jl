# Test with linear chirp

rate = 16000.
duration = 2
samples = round(Int,duration*rate)
y = [ t <= 8*samples/14 ? sin(1000*2*π*t/rate+1000*t/rate*2*π*t/rate) : 0. for t in 1:samples ]

mkpath("linear-chirp-results")
cd("linear-chirp-results")
wavwrite(y, "linear_chirp.wav", Fs = rate)

@time SS = stft(y, 500, 450, fs = rate)

@time Lm = lift(SS, threshold = 10, N=30)

KK=20
AA=110
BB=50
CC=250


τ = KK*step(time(Lm))
@time k = kernel_computation(normalize(freq(Lm)), slopes(Lm), τ, n = 20);

@time W = wc_delay(Lm, AA, BB, CC, K=k)

save_result(W, AA, BB, CC, KK)

