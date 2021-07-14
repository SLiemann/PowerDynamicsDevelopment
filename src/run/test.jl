using Plots
#using PowerDynamics
#using DifferentialEquations
#using ModelingToolkit
using FFTW
using ModelingToolkit



f(x;freq=1) = sin(x*2*pi*freq)
abtast = 1e-3
t =collect(0:abtast:5)
signal = f.(t).+f.(t,freq=5) .+0.5*f.(t,freq=3)

function DFT(signal,t)
    N = length(t)
    Ts = t[end]/N  #it is assumed that measured point are equally distributed
    # Fourier Transform of it
    F = fft(signal) |> fftshift
    freqs = fftfreq(N, 1.0/Ts) |> fftshift
    return F, freqs
end
tmp_F, tmp_freqs = DFT(signal,t)
plot(tmp_freqs, abs.(tmp_F), title = "Spectrum", xlim=(0, +10))

@variables x, y, z
eqs = [y ~ atan(x),
       z ~ sin(x) + cos(x)]
expand_derivatives.(Differential(x).(eqs))
expand_derivatives(Differential(x)(sin(x)+cos(x)))
xlims!((0.75,1))
# Number of points
N = 2^14
# Sample period
Ts = 1 / (1.1 * N)
# Start time
t0 = 0
tmax = t0 + N * Ts
# time coordinate
t = t0:Ts:tmax

# signal
signal = sin.(2π * 60 .* t) # sin (2π f t)

# Fourier Transform of it
F = fft(signal) |> fftshift
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

# plots
time_domain = plot(t, signal, title = "Signal")
freq_domain = plot(freqs, abs.(F), title = "Spectrum", xlim=(-100, +100))
plot(time_domain, freq_domain, layout = 2)
