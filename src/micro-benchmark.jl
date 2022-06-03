using BenchmarkTools
function cumuSum(n)
    A = rand(Float64, n)
    acc = 0
    for i in A
        acc += i
    end
end

n = 100
while n<=100000000
    @btime cumuSum(n)
    global n*=10
end

