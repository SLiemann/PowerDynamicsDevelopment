using ModelingToolkit
using DifferentialEquations
using PowerDynamics
using FFTW

function getPreFaultVoltages(pg::PowerGrid,ic_prefault::Array{Float64,1},ic_endfault::Array{Float64,1})
    ind = getVoltageSymbolPositions(pg)
    ic_endfault[ind] = ic_prefault[ind]
    return ic_endfault
end

function getPreFaultAlgebraicStates(pg::PowerGrid,ic_prefault::Array{Float64,1},ic_endfault::Array{Float64,1})
    prob = ODEFunction(rhs(pg))
    ind = findall(x-> iszero(x),diag(prob.mass_matrix))
    ic_endfault[ind] = ic_prefault[ind]
    return ic_endfault
end

function getVoltageSymbolPositions(pg::PowerGrid)
    n = map(x->PowerDynamics.variable_index(pg.nodes,x,:u_r),collect(keys(pg.nodes)))
    append!(n,map(x->PowerDynamics.variable_index(pg.nodes,x,:u_i),collect(keys(pg.nodes))))
end

function getSymbolPosition(pg::PowerGrid,syms::Array{Symbol,1})
    return sort(indexin(syms,rhs(pg).syms))
end

function AddNaNsIntoSolution(pg1::PowerGrid,pg2::PowerGrid,sol)
    key1 = collect(keys(pg1.nodes))
    key2 = collect(keys(pg2.nodes))
    n = findfirst(x->x==nothing,indexin(key1,key2))
    if ~isnothing(n)
        len_bus = length(symbolsof(pg1.nodes[key1[n]]))
        first_ind = PowerDynamics.variable_index(pg1.nodes,key1[n],:u_r)
        start_u = findfirst(x-> x== length(sol.u[end]),length.(sol.u))
        for i in start_u:length(sol.u)
            splice!(sol.u[i],first_ind:first_ind-1,NaN*zeros(len_bus))
            for j in sol.k[i]
                splice!(j,first_ind:first_ind-1,NaN*zeros(len_bus))
            end
        end
    end
    return sol
end

function PiModel(y::Complex{Num}, y_shunt_km, y_shunt_mk, t_km, t_mk)
    Π = Array{Union{Complex{Float64},Complex{Num}}}(undef,2,2)
    Π[1, 1] = - abs2(t_km) * (y + y_shunt_km) # Our sign convention is opposite for the source of the edge
    Π[1, 2] = conj(t_km) * t_mk * y # Our sign convention is opposite for the source of the edge
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
end

function PiModel(y, y_shunt_km, y_shunt_mk, t_km, t_mk::Num)
    Π = Array{Union{Complex{Float64},Complex{Num}}}(undef,2,2)
    Π[1, 1] = - abs2(t_km) * (y + y_shunt_km) # Our sign convention is opposite for the source of the edge
    Π[1, 2] = conj(t_km) * t_mk * y # Our sign convention is opposite for the source of the edge
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
end

function PiModel(y, y_shunt_km, y_shunt_mk, t_km::Num, t_mk)
    Π = Array{Union{Complex{Float64},Complex{Num}}}(undef,2,2)
    Π[1, 1] = - abs2(t_km) * (y + y_shunt_km) # Our sign convention is opposite for the source of the edge
    Π[1, 2] = conj(t_km) * t_mk * y # Our sign convention is opposite for the source of the edge
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
end

function CalcEigenValues(pg::PowerGrid, p::Array{Float64,1}; output::Bool = false, plot::Bool = false)
  mtsys = GetMTKSystem(pg, (0.0, 1.0), p)
  Fx, Fy, Gx, Gy = GetSymbolicFactorizedJacobian(mtsys)
  Fxf, Fyf, Gxf, Gyf = [
    Substitute(f, mtsys.defaults) for
    f in [Fx, Fy, Gx, Gy]
  ]
  Af = Fxf - Fyf * inv(Gyf) * Gxf
  EW = eigvals(Af)

  x,y = GetSymbolicStates(mtsys)
  index = indexin(x, states(mtsys))
  syms = rhs(pg).syms[index]
  if output
    display("|ID | Real-part | Imag-part | Frequency | Damping Time Constant |")
    for (ind, ew) in enumerate(EW)
      display(
        "| $(syms[ind])) | $(round(real(ew),digits =3)) | $(round(imag(ew),digits = 3)) | $(round(abs(imag(ew))/2/pi,digits =3)) | $(round(1.0/abs(real(ew)),digits =3)) |",
      )
    end
  end
  if plot
      scatter([real(EW)[1]],[imag(EW)[1]], legend = :outertopright,label = String(syms[1]))
      for i in 2:length(EW)
          scatter!([real(EW)[i]],[imag(EW)[i]],label = String(syms[i]))
      end
      ylims!((-maximum(imag.(EW))*1.1,maximum(imag.(EW))*1.1))
      min_real_ew = minimum(real.(EW))
      plot!([min_real_ew;0.0],[min_real_ew*20.0;0.0], linestyle=:dash,linecolor = :red, label = "5 % damping")
      plot!([min_real_ew;0.0],[-min_real_ew*20.0;0.0], linestyle=:dash,linecolor = :red, label = nothing)
      plot!([min_real_ew;0.0],[min_real_ew*10.0;0.0], linestyle=:dash,linecolor = :blue, label = "10 % damping")
      plot!([min_real_ew;0.0],[-min_real_ew*10.0;0.0], linestyle=:dash,linecolor = :blue, label = nothing)
      plot!([min_real_ew;0.0],[min_real_ew*5.0;0.0], linestyle=:dash,linecolor = :green, label = "20 % damping")
      display(plot!([min_real_ew;0.0],[-min_real_ew*5.0;0.0], linestyle=:dash,linecolor = :green, label = nothing))
  end
  return EW
end

function DFT(signal,t)
    N = length(t)
    Ts = t[end]/N  #it is assumed that measured point are equally distributed
    # Fourier Transform of it
    F = fft(signal) |> fftshift
    freqs = fftfreq(N, 1.0/Ts) |> fftshift
    return F, freqs
end

function ExtractResult(pgsol::PowerGridSolution, sym::Symbol)
    sol = pgsol.dqsol
    ind = indexin([sym],rhs(pgsol.powergrid).syms)[1]
    return pgsol.dqsol[ind,:]
end

function ExtractResult(pgsol::PowerGridSolution, bus::String, sym::Symbol)
    keys_node = collect(keys(pgsol.powergrid.nodes))
    ind = indexin([bus],keys_node)[1]
    if isnothing(ind)
        error("$bus is not a valid bus key")
    else
        sym = Symbol(string(sym,"_",ind))
        return ExtractResult(pgsol,sym)
    end
end
