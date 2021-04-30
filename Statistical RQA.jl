
using RecurrenceAnalysis, Statistics



function LAM_Probabilistic(mapp::Array{Float64, 1}; n=n, q=2, lmin=2, samples=Int(n/10), e=0.3)
    #=
        This function will calculate the determinism using a probabilistical approach. It will calculate various DET and average over them.
        Input:
            mapp: time series, data.
            n: the size of the time series
            q: the microstate that should be used
            lmin: the minimum line length you want to consider
            samples: the sampling rule you want to test (the suggested ]fast sampling rule is: floot(n/q) - note this is the default)
            e: the threshold parameter you want
        Outpub:
            DET: the value of the determinism
    =#

    LAM = zeros(Float64, samples)

    for s in 1:samples
        s = Int(s)
        a, b = rand(1:n-q+1), rand(1:n-q+1)
        micro = CrossRecurrenceMatrix(mapp[a:a+q-1], mapp[b:b+q-1], e)
        LAM[s] = laminarity(micro, theiler=0)
    end

    return mean(LAM)
end

function TT_Probabilistic(mapp::Array{Float64, 1}; n=n, q=2, lmin=2, samples=Int(n/10), e=0.3)
    #=
        This function will calculate the determinism using a probabilistical approach. It will calculate various DET and average over them.
        Input:
            mapp: time series, data.
            n: the size of the time series
            q: the microstate that should be used
            lmin: the minimum line length you want to consider
            samples: the sampling rule you want to test (the suggested ]fast sampling rule is: floot(n/q) - note this is the default)
            e: the threshold parameter you want
        Outpub:
            DET: the value of the determinism
    =#

    TT = zeros(Float64, samples)

    for s in 1:samples
        s = Int(s)
        a, b = rand(1:n-q+1), rand(1:n-q+1)
        micro = CrossRecurrenceMatrix(mapp[a:a+q-1], mapp[b:b+q-1], e)
        TT[s] = trappingtime(micro, theiler=0)
    end

    return mean(TT)
end

function DET_Probabilistic(mapp::Array{Float64, 1}; n=n, q=2, lmin=2, samples=0, e=0.3)
    #=
        This function will calculate the determinism using a probabilistical approach. It will calculate various DET and average over them.
        Input:
            mapp: time series, data.
            n: the size of the time series
            q: the microstate that should be used
            lmin: the minimum line length you want to consider
            samples: the sampling rule you want to test (the suggested ]fast sampling rule is: floot(n/q) - note this is the default)
            e: the threshold parameter you want
        Outpub:
            DET: the value of the determinism
    =#

    if samples == 0
        samples = Int(floor(n/q))
    end

    DET = zeros(Float64, Int(floor(samples)))

    for s in 1:samples
        s = Int(s)
        a, b = rand(1:n-q+1), rand(1:n-q+1)
        micro = CrossRecurrenceMatrix(mapp[a:a+q-1], mapp[b:b+q-1], e)
        DET[s] = determinism(micro)
    end

    return mean(DET)
end

function L_Probabilistic(mapp::Array{Float64, 1}; n=n, q=2, lmin=2, samples=0, e=0.3)
    #=
        This function will calculate the determinism using a probabilistical approach. It will calculate various DET and average over them.
        Input:
            mapp: time series, data.
            n: the size of the time series
            q: the microstate that should be used
            lmin: the minimum line length you want to consider
            samples: the sampling rule you want to test (the suggested ]fast sampling rule is: floot(n/q) - note this is the default)
            e: the threshold parameter you want
        Outpub:
            DET: the value of the determinism
    =#

    if samples == 0
        samples = Int(floor(n/q))
    end

    L = zeros(Float64, Int(floor(samples)))

    for s in 1:samples
        s = Int(s)
        a, b = rand(1:n-q+1), rand(1:n-q+1)
        micro = CrossRecurrenceMatrix(mapp[a:a+q-1], mapp[b:b+q-1], e)
        L[s] = average_l(micro)
    end

    return mean(L)
end
