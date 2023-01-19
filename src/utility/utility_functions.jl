using ModelingToolkit
using DifferentialEquations
using PowerDynamics
using FFTW
using LinearAlgebra
import PlotlyJS: plot

function getStateIndex(pg::PowerGrid,node_str::String,sym::Symbol)
    ind = 0
    node = pg.nodes[node_str]
    for (ind_tmp,val) in enumerate(pg.nodes)
        if val[2] == node
            ind = ind + indexin([sym],symbolsof(node))[1]
            break
        end
       ind += length(symbolsof(val[2]))
    end
    return ind
end

function getPreFaultVoltages(pg::PowerGrid,ic_prefault::Array{Float64,1},ic_endfault::Array{Float64,1})
    ind = getVoltageSymbolPositions(pg)
    ic_endfault[ind] = ic_prefault[ind]
    return ic_endfault
end

function getPreFaultVoltages(pg::PowerGrid,ic_prefault,ic_endfault)
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

function getComplexBusVoltage(pg::PowerGrid,ic::Array{Float64,1})
        ind = sort(getVoltageSymbolPositions(pg))
        Uc = Complex.(zeros(Int(length(ind)/2)))
        for (i,val) in enumerate(ind[1:2:length(ind)])
            Uc[i] = ic[val] + 1im*ic[val+1]
        end
        return Uc
end

function getSymbolPosition(pg::PowerGrid,syms::Symbol)
    ind =  indexin([syms],rhs(pg).syms)[1]
    if isnothing(ind)
        error(syms," is not a valid symbol")
    end
    return ind
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

# The following PiModel function is copied from the actual PowerDynamics package, Michael
function PiModel(y, y_shunt_km, y_shunt_mk, t_km, t_mk)
    Π = Matrix{Any}(undef,2,2)#zeros(Complex{Float64}, 2, 2)
    Π[1, 1] = - abs2(t_km) * (y + y_shunt_km) # Our sign convention is opposite for the source of the edge
    Π[1, 2] = conj(t_km) * t_mk * y # Our sign convention is opposite for the source of the edge
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
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

function CalcEigenValues(pg::PowerGrid, p; output::Bool = false, plot::Bool = false)
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
      display(plot!([min_real_ew;0.0],[-min_real_ew*5.0;0.0], linestyle=:dash,linecolor = :green, label = nothing));
  end
  return EW
end

function DFTplot(signal,t,fmax;norm=true)
    F,freqs = DFT(signal,t,norm_to_fundamental=norm)
    display(bar(freqs,abs.(F),xlims=(0,fmax)))
end

function DFT(signal,t;norm_to_fundamental=false)
    N = length(t)
    Ts = t[end]/N  #it is assumed that measured point are equally distributed
    # Fourier Transform of it
    F = fft(signal) |> fftshift
    freqs = fftfreq(N, 1.0/Ts) |> fftshift
    if norm_to_fundamental
        ind = indexin([0.0],freqs)[1]+1
        F ./= abs.(F[ind])
    end
    return F, freqs
end

function RMS(signal,time;fund = 50)
    RMS = zeros(length(signal))
    dt = time[end]/length(time) # fixed step-size assumed
    L = Int(floor(1/fund/dt))
    zaehler = Int(floor(1/fund/dt))
    ind = 0;
    for i=0:1:length(signal)-L
        s = signal[i+1:i+L]
        t = time[i+1:i+L]
        t = t.-t[1]
        val, freq = DFT(s,t)
        ind = i==0 ? findmin(abs.(freq.-fund))[2] : ind
        RMS[zaehler] = abs(val[ind]);
        zaehler = zaehler + 1;
    end
    RMS./(sqrt(2)*L/2)
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

function SavePlotApprTrajectories(pg_tmp,sol_or,sol_per,sens,par_num,or_value,new_value,labels_p; state_sym = :i_abs)
    sol_appr = deepcopy(sol_or.dqsol)
    for (ind,val) in enumerate(collect(eachcol(sens[par_num])))
        sol_appr.u[ind+1] .+= val*(new_value-or_value)
    end
    pgsol_tmp = PowerGridSolution(sol_appr,pg_tmp)
    title_str  = labels_p[par_num] * ": from " * string(or_value) * " to " * string(new_value)
    plot(sol_or,"bus3",state_sym, label = "Original - " * string(state_sym), title = title_str)
    plot!(sol_per,"bus3",state_sym, label = "Real perturbed - " * string(state_sym))
    display(plot!(pgsol_tmp,"bus3",state_sym, label = "Approximated - " * string(state_sym),linestyle = :dash))
    savefig("C:\\Users\\liemann\\Desktop\\Julia_figures\\"*labels_p[par_num] *".svg")
end

function PlotApproTrajectories(
        pg::PowerGrid,
        sol_or::PowerGridSolution,
        sol_per::PowerGridSolution,
        sensis::Vector{Array{Float64}},
        par_num::Int64,
        or_value::Float64,
        new_value::Float64,
        labels_p::Vector{String},
        state::Symbol;
        bus = "bus4"
    )
    sol_appr = CalcApprSolution(sol_or,sensis,par_num,or_value,new_value)
    title_str  = labels_p[par_num] * ": from " * string(or_value) * " to " * string(new_value)
    plot(sol_or,bus,state, label = "Original - " * string(state), title = title_str)
    plot!(sol_per,bus,state, label = "Real perturbed - " * string(state))
    display(plot!(sol_appr,bus,state, label = "Approximated - " * string(state),linestyle = :dash))
end

function CalcApprSolution(sol_or::PowerGridSolution,sensis::Vector{Array{Float64}},par_num::Int64,or_value::Float64,new_value::Float64)
    sol_appr = deepcopy(sol_or)
    for (ind,val) in enumerate(collect(eachcol(sensis[par_num])))
        sol_appr.dqsol.u[ind+1] .+= val*(new_value-or_value)
    end
    return sol_appr
end

function SaveComparingTrajectoryPlots(pg_tmp,ic_or,pgsol_or,sens,labels_p,delta)
    params_tmp = GFC_params()
    for (ind,val) in enumerate(params_tmp)
        params_new = deepcopy(params_tmp)
        params_new[ind] += delta*params_new[ind]
        prob_tmp = ODEProblem(rhs(pg_tmp),ic_or,(0.0,3.0),params_new)
        try
            pgsol_per,evr = simGFC(prob_tmp)
            SavePlotApprTrajectories(pg_tmp,pgsol_or,pgsol_per,sens,ind,params[ind],params_new[ind],labels_p)
        catch
            @warn "Could not calc: " *labels_p[ind]
        end
    end
end

function GetVoltageSensis(sensis,ind_ur_state,ind_ui_state)
    sensi_ur = Array{Float64}(undef,size(sensis[1])[2],length(sensis));
    sensi_ui = Array{Float64}(undef,size(sensis[1])[2],length(sensis));

    for (ind,val) in enumerate(sensis)
        sensi_ur[:,ind] = val[ind_ur_state,:]
        sensi_ui[:,ind] = val[ind_ui_state,:]
    end
    return sensi_ur, sensi_ui
end

function GetAbsVoltageSensis(pg_sol::PowerGridSolution,ur::Symbol,ui::Symbol,sensis::Vector{Array{Float64}},param::Vector{Int64},Δp::Vector{Float64})
    #Idea is: |U_appr| = |U_or| + Z   with |U_or| = sqrt(Ur_or^2 + Ui_or^2)
    #         U_r_appr = Ur_or  + xx0*Δp (same for U_i_appr)
    #         |U_appr| = sqrt(U_r_appr^2 + U_r_appr^2)
    #next    |U_appr|^2 = (|U_or| + Z)^2 solve for Z

    ur_or = ExtractResult(pg_sol,ur)[1:end-1]
    ui_or = ExtractResult(pg_sol,ui)[1:end-1]

    ind = indexin([ur,ui],rhs(pg_sol.powergrid).syms)
    sens_ur,sens_ui = GetVoltageSensis(toll_tap,ind[1],ind[2])
    z1 = Array{Float64}(undef,length(pg_sol.dqsol.t[1:end-1]),length(param))
    z2 = Array{Float64}(undef,length(pg_sol.dqsol.t[1:end-1]),length(param))
    for (ind,val) in enumerate(param)
        z1[:,ind],z2[:,ind] = CalcAbsVoltageSens(ur_or,ui_or,sens_ur,sens_ui,val,Δp[ind])
    end
    return z1,z2
end

function CalcAbsVoltageSens(ur_or,ui_or,sens_ur,sens_ui,param,Δp)
    s_ur = sens_ur[:,param]*Δp
    s_ui = sens_ur[:,param]*Δp

    p = 2*sqrt.(ur_or.^2.0 .+ ui_or.^2.0)
    q = -(2 .*(ur_or.*s_ur + ui_or.*s_ui) .+ s_ur.^2 .+ s_ui.^2)

    z1 = -p./2 + sqrt.((p./2).^2 .- q)
    z2 = -p./2 - sqrt.((p./2).^2 .- q)
    return z1,z2
end

function CalcAbsVoltageSens_v2(ur_or,ui_or,sens_ur,sens_ui;Δp=1.0)
    s_ur = sens_ur*Δp
    s_ui = sens_ur*Δp

    p = 2*sqrt.(ur_or.^2.0 .+ ui_or.^2.0)
    q = -(2 .*(ur_or.*s_ur + ui_or.*s_ui) .+ s_ur.^2 .+ s_ui.^2)

    z1 = -p./2 + sqrt.((p./2).^2 .- q)
    z2 = -p./2 - sqrt.((p./2).^2 .- q)
    return z1,z2
end

function getallParameters(node)
    cond = true
    i = 1
    param = Vector{Float64}()
    while cond
        try
            tmp = getfield(node,i)
            append!(param,tmp)
            i += 1
        catch
            cond = false
        end
    end
    param[1:end-1] #last entry is Y_n from PowerDynamics
end

function Sol2DF(pgsol::PowerGridSolution)
    sol = DataFrame(pgsol.dqsol)

    for (ind,val) in enumerate(pgsol.powergrid.nodes)
        ur = ExtractResult(pgsol,Symbol("u_r_"*string(ind)))
        ui = ExtractResult(pgsol,Symbol("u_i_"*string(ind)))
        u = sqrt.(ur.^2 + ui.^2)
        sol[!,Symbol("uabs_"*string(ind))] = u;
    end
    return sol
end


function plotallvoltages(pgsol::PowerGridSolution)
    newdf = Sol2DF(pgsol)
    p = Vector{GenericTrace}()
    for (ind,val) in enumerate(pgsol.powergrid.nodes)
        tmp = scatter(x=newdf.timestamp,y=newdf[:,"uabs_"*string(ind)],name=val[1])
        push!(p,tmp)
    end
    display(plot(p))
    return p
end

function myplot(pgsol::PowerGridSolution,sym::Symbol;y_norm=1.0,y_bias = 0.0)
    ind = findfirst(x->x==sym,collect(rhs(pgsol.powergrid).syms))
    t = pgsol.dqsol.t
    y =  pgsol.dqsol[ind,:]./y_norm .+ y_bias
    sc = scatter(x=t,y=y,name=String(sym))
    #display(plot(sc))
    return sc
end

function myplot(pgsol::PowerGridSolution,bus::String,sym::Symbol;y_norm=1.0,y_bias = 0.0)
    ind_bus = findfirst(x->x==bus,collect(keys(pgsol.powergrid.nodes)))
    sym = Symbol(string(sym)*"_"*string(ind_bus))
    myplot(pgsol,sym,y_norm=y_norm,y_bias=y_bias)
end

function myplot(pgsol::Vector{PowerGridSolution},bus::String,sym::Symbol;y_norm=1.0,y_bias = 0.0)
    p = Vector{GenericTrace}()
    for i in pgsol
        tmp  =myplot(i,bus,sym,y_norm=y_norm,y_bias=y_bias)
        push!(p,tmp)
    end
    #display(plot(p))
    return p
end

function plotv(pgsol,bus::String)
    ind = findfirst(x->x==bus,collect(keys(pgsol.powergrid.nodes)))
    ur = ExtractResult(pgsol,Symbol("u_r_"*string(ind)))
    ui = ExtractResult(pgsol,Symbol("u_i_"*string(ind)))
    t = pgsol.dqsol.t
    y =  sqrt.(ur.^2 + ui.^2)
    sc = scatter(x=t,y=y,name=bus)
    return sc
end

function plotv(pgsol,bus::Vector{String})
    p = Vector{GenericTrace}()
    for (ind,val) in enumerate(bus)
        tmp = plotv(pgsol,val)
        push!(p,tmp)
    end
    return p
end

function plotv(pgsol::Vector{PowerGridSolution},bus::String)
    p = Vector{GenericTrace}()
    for i in pgsol
        tmp  = plotv(i,bus)
        push!(p,tmp)
    end
    return p
end

function DetermineBoundary(XR,Rv,Xv)
    r = Vector{Float64}()
    x = Vector{Float64}()
    flag_dir = "left"
    flag_end = false
    ind = [(-1,0);(0,1);(1,0);(0,-1)]

    XR = reverse(XR)'
    Rv = reverse(Rv)
    Xv = reverse(Xv)

    lenx = size(XR)[1]
    lenr = size(XR)[2]

    #starting indices
    indr = 1
    indx = findfirst(k->k==1,XR[:,1])

    #saving first entry
    push!(r,Rv[indr])
    push!(x,Xv[indx])

    function next_ind(start_ind,r,x)
        for i=start_ind:start_ind+3
            if mod(i,4) == 0
                dx = ind[4][1]
                dr = ind[4][2]
            else
                dx = ind[mod(i,4)][1]
                dr = ind[mod(i,4)][2]
            end

            if indr+dr == 0 || indr+dr > lenr || indx+dx > lenx 
                #skip boundarys
            elseif indx + dx == 1  && XR[indx,indr] == 1 
                indr = indr + dr
                indx = indx + dx
                push!(r,Rv[indr])
                push!(x,Xv[indx])
                indr = indr -1
                
                while indr != 0 && XR[indx,indr] == 1  # go x-axis to zero 
                    push!(r,Rv[indr])
                    push!(x,Xv[indx])
                    indr = indr -1
                end
                flag_end = true
                break;
            elseif XR[indx+dx,indr+dr] != 1 
                #skip unstable case
            elseif XR[indx + dx,indr + dr] == 1
                indr = indr + dr
                indx = indx + dx
                push!(r,Rv[indr])
                push!(x,Xv[indx])

                if mod(i,4) == 1
                   flag_dir = "bottom"
                   break;
                elseif mod(i,4) == 2
                   flag_dir = "left"
                   break;
                elseif mod(i,4) == 3
                    flag_dir = "top"
                    break;
                elseif mod(i,4) == 0
                    flag_dir = "right"
                    break;
                else
                    error("error2")
                end
            end 
        end
        return flag_end,flag_dir,r,x
    end    
    counter = 0;

    while flag_end == false && counter < 1e4
        if flag_dir == "top"
            flag_end,flag_dir,r,x = next_ind(2,r,x)
        elseif flag_dir == "left"
            flag_end,flag_dir,r,x = next_ind(1,r,x)
        elseif flag_dir == "bottom"
            flag_end,flag_dir,r,x = next_ind(4,r,x)
        elseif flag_dir == "right"
            flag_end,flag_dir,r,x = next_ind(3,r,x)
        else
            error("error")
        end
        counter = counter + 1
    end
    if counter == 1e4
        error("Counter reached 1e4")
    end
    return r,x
end

function CompareXRResults(obj::Array{String})
    sc = Vector{GenericTrace}()
    for i in obj
        data = load(i)
        Rv= data["Rverlauf"];
        Xv = data["Xverlauf"];
        XRv= data["XR"];
        x,y = DetermineBoundary(XRv,Rv,Xv)
        tmp_sc = scatter(x=x,y=y,name=i)
        push!(sc,tmp_sc)
    end
    display(plot(sc))
end