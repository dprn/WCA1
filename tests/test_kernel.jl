function test_threshold(f, s; p = KernParams(.1,1,1e-3))
#     println("Testing ", tol(p))
    k = Kern(f,s, p)
    m = 0
    for i in 1:length(f), j in 1:length(s)
        kt = [ WCA1.K(p, (f[i], s[j]), (f[r1], s[r2])) for r1 in 1:length(f), r2 in 1:length(s) ]
        for I in k[i,j]
            kt[I[1]] = 0
        end
        
        if any( kt .> tol(p) )
            println("Failed at", i, " ", j)
            return false
        end
    end
    
    # println("Everything ok (tol = ", tol(p), ").")
    return true
end

f = -5:.1:5
s = -5:.1:5
τ = .1
b = 1

@testset "Kernel threshold test" begin
    @test test_threshold(f, s, p = KernParams(τ, b, 1e-3))
    @test test_threshold(f, s, p = KernParams(τ, b, 1e-4))
    @test test_threshold(f, s, p = KernParams(τ, b, 1e-5))
    @test test_threshold(f, s, p = KernParams(τ, b, 1e-6))
end

#=test_threshold(f, s, p = KernParams(τ, b, 1e-3))
test_threshold(f, s, p = KernParams(τ, b, 1e-4))
test_threshold(f, s, p = KernParams(τ, b, 1e-5))
test_threshold(f, s, p = KernParams(τ, b, 1e-6))=#