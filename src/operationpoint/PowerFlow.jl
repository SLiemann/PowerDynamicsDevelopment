using PowerDynamics
using PowerDynamics: guess
using LightGraphs #incidence_matrix
using Roots #for Iwamoto multiplier
using LinearAlgebra

#Pi models for nodal admittance matrice
PiModel(L::PiModelLine) = PiModel(L.y,L.y_shunt_km,L.y_shunt_mk,1,1)
PiModel(T::Transformer) = PiModel(T.y,0,0,T.t_ratio,1)
PiModel(S::StaticLine)  = PiModel(S.Y,0,0,1,1)
PiModel(R::RLLine)      = PiModel(1/(R.R+R.ω0*R.L),0,0,1,1)
function PiModel(T::StaticPowerTransformer)
    u_LV = T.U_LV
    u_HV = T.U_HV
    if T.tap_side == "LV"
        u_LV = (1+T.tap_pos*T.tap_inc/100.)*T.U_LV
    elseif T.tap_side == "HV"
        u_HV = (1+T.tap_pos*T.tap_inc/100.)*T.U_HV
    else
        error("Can not interprete tap_side (HV/LV): $tap_side")
    end
    ü     = u_HV/u_LV
    I_r   = T.S_r/(u_HV*sqrt(3)) #rated current

    #Calculatiing leakage reactance Xa and winding resistance Ra
    Ra = 0.
    Xa = 0.
    if T.XR_ratio == 0
        Ra = T.uk*u_HV/(sqrt(3)*I_r);
        Xa = 0.
    elseif T.XR_ratio == Inf
        Ra = 0.
        Xa = T.uk*u_HV/(sqrt(3)*I_r);
    else
        Ra = sqrt(((T.uk*u_HV/(I_r*sqrt(3)))^2)/(1+T.XR_ratio^2));
        Xa = T.XR_ratio*Ra;
    end
    Ya = 1/(0.5*(Ra+1im*Xa)) #it is assumed that losses are equally (0.5) distributed over both sides
    Ybs = Ya; #Ybs should be changed, if losses are not equally distributed

    #Calculating magnetising reactance Xm and core resistance Rfe from iron losess
    Zm  = u_HV/(sqrt(3)*T.i0/100.0*I_r) #no-load currents depends on complete magnetising impedance
    Rfe = u_HV*u_HV/(T.Pv0);
    Xm = Inf
    if T.i0 != 0
        Xm  = 1/(sqrt(1/Zm^2 - 1/Rfe^2))
    end
    Ym  = 1/Rfe + 1/(1im*Xm)

    Ybase = 1/(T.Ubase^2/T.Sbase) #this could be changed, if global base values are available
    Y     = 1/(Ya+Ybs+Ym)./Ybase
    Y = Y.*PiModel(Ya*Ybs, Ya*Ym, Ybs*Ym,1,ü)
    return Y
end

# NodeTypes: 0 = Slack, 1 = PV, 2 = PQ
NodeType(S::SlackAlgebraic) = 0
NodeType(F::FourthOrderEq)  = 1
NodeType(F::FourthOrderEqExciterIEEEDC1A)  = 1
NodeType(F::FourthOrderEqGovernorExciterAVR)  = 1
NodeType(F::FourthOrderEqGovernorIEEEG1)  = 1
NodeType(S::SwingEq)  = 1
NodeType(L::PVAlgebraic) = 1
NodeType(V::VSIMinimal) = 2
NodeType(V::VSIVoltagePT1) = 2
NodeType(L::PQAlgebraic) = 2
NodeType(L::VoltageDependentLoad) = 2
NodeType(L::ExponentialRecoveryLoad)  = 2
NodeType(L::CSIMinimal)  = 2

#note: only loads are treated with voltage depency and are called every iteration
PowerNodeLoad(S::SlackAlgebraic,U) = 0. #treated as generation
PowerNodeLoad(F::FourthOrderEq,U) = 0. #treated as generation
PowerNodeLoad(F::FourthOrderEqExciterIEEEDC1A,U)  = 0. #treated as generation
PowerNodeLoad(F::FourthOrderEqGovernorExciterAVR,U)  = 0. #treated as generation
PowerNodeLoad(F::FourthOrderEqGovernorIEEEG1,U)  = 0. #treated as generation
PowerNodeLoad(S::SwingEq,U)  = 0. #treated as generation
PowerNodeLoad(V::VSIMinimal,U) = complex(V.P,(abs(U)-V.V_r)/V.K_Q+V.Q)   #treated as load with changed sign to include voltage dependency
PowerNodeLoad(V::VSIVoltagePT1,U) = complex(V.P,(abs(U)-V.V_r)/V.K_Q+V.Q)   #treated as load with changed sign to include voltage dependency
PowerNodeLoad(L::PVAlgebraic,U) = -L.P  #treated as load
PowerNodeLoad(L::PQAlgebraic,U) = -complex(L.P,L.Q) #treated as load
PowerNodeLoad(L::VoltageDependentLoad,U) = -complex(L.P, L.Q) * (L.A * abs(U)^2 + L.B * abs(U) + 1 - L.A - L.B)
PowerNodeLoad(L::ExponentialRecoveryLoad,U)  = -(L.P0*((abs(U)/L.V0)^L.Nps) + 1im*L.Q0*((abs(U)/L.V0)^L.Nqs))
PowerNodeLoad(L::CSIMinimal,U)  = -U*conj(L.I_r)

#generation is voltage independent, otherwise it has to be called every iteration
PowerNodeGeneration(S::SlackAlgebraic) = 0.
PowerNodeGeneration(F::FourthOrderEq) = F.P
PowerNodeGeneration(F::FourthOrderEqExciterIEEEDC1A)  = F.P
PowerNodeGeneration(F::FourthOrderEqGovernorExciterAVR)  = F.P
PowerNodeGeneration(F::FourthOrderEqGovernorIEEEG1)  =  F.P
PowerNodeGeneration(S::SwingEq)  = S.P
PowerNodeGeneration(V::VSIMinimal) = 0. #treated as load with changed sign to include voltage dependency
PowerNodeGeneration(V::VSIVoltagePT1) = 0. #treated as load with changed sign to include voltage dependency
PowerNodeGeneration(L::PVAlgebraic) = 0. #treated as load
PowerNodeGeneration(L::PQAlgebraic) = 0. #treated as load
PowerNodeGeneration(L::VoltageDependentLoad) = 0. #treated as load
PowerNodeGeneration(L::ExponentialRecoveryLoad)  = 0. #treated as load
PowerNodeGeneration(L::CSIMinimal)  = 0. #treated as load


function PowerFlowClassic(pg::PowerGrid, U_r_nodes::Vector{Float64}; Ubase::Float64=380e3, ind_sl::Int64 = 0,max_tol::Float64 = 1e-5,iter_max::Int64  = 30,iwamoto::Bool =false, Qmax = -1, Qmin = -1)
    number_nodes = length(pg.nodes); #convenience
    nodetypes = NodeType.(values(pg.nodes))
    if !isempty(findall(x-> x==0, nodetypes)) #if there is no SlackAlgebraic
        ind_sl = findall(x-> x==0, nodetypes)[1] #set passed value
        @info "Reference node is bus no. $ind_sl"
    end
    ind_PV_or = findall(x-> x==1, nodetypes)
    ind_PQ_or = findall(x-> x==2, nodetypes)

    ind_PQ = ind_PQ_or #copy is needed for checking Qlimits
    ind_PV = ind_PV_or #copy is needed for checking Qlimits

    #to calculate Ykk a vector with all SI voltages of the nodes is needed (U_r_nodes)
    Ykk = NodalAdmittanceMatrice(pg,U_r_nodes,Ubase);

    U = ones(number_nodes,1);
    δ = CalcδStartValues(pg,Ykk,ind_sl);

    Ykk_abs = abs.(Ykk);
    θ   = angle.(Ykk);
    Pn  = similar(U);
    Qn  = similar(U);
    S_node_gen = Array{Complex{Float64},2}(undef,number_nodes,1)
    S_node_load = Array{Complex{Float64},2}(undef,number_nodes,1)
    #PowerNodeGeneration is called only once, since it should only change,
    #when power limits are reached or voltage dependency should be included
    for (ind,val) in enumerate(values(pg.nodes)) S_node_gen[ind] = PowerNodeGeneration(val) end

    #start of the iterationen
    for iter in 1:iter_max+1
        for i in 1:number_nodes #calculate node powers with current voltage
            Pn[i] = sum(U[i].*U.*Ykk_abs[:,i].*cos.(δ[i].-δ.-θ[:,i]))
            Qn[i] = sum(U[i].*U.*Ykk_abs[:,i].*sin.(δ[i].-δ.-θ[:,i]))
        end
        #update load powers
        for (ind,val) in enumerate(values(pg.nodes)) S_node_load[ind] = PowerNodeLoad(val,U[ind]*exp(1im*δ[ind])) end

        S_node = S_node_gen + S_node_load #total power at a node

        ΔP = real(S_node) - Pn
        ΔQ = imag(S_node) - Qn

        ΔP = ΔP[1:end .!= ind_sl, :]; #delete slack
        ΔQ = ΔQ[setdiff(1:end, [ind_sl; ind_PV]), :]; #delete slack and PV nodes

        error_ = norm([ΔP;ΔQ],Inf)
        if error_ < max_tol
            @info "Power flow converged in $iter iterations"
            ic_guess = guess.(values(pg.nodes),U.*exp.(1im*δ))
            return U,δ*180/pi,vcat(ic_guess...) # power flow converged
        elseif iter == iter_max
             @warn "Power flow reached max. iteration ($iter) and has not converged."
             ic_guess = guess.(values(pg.nodes),U.*exp.(1im*δ))
             return vcat(ic_guess...) #max iteration reached
        end

        #get load flow Jacobian
        J = CalculatePolarLoadFlowJacobian(U,δ,Ykk)
        del_sl = [ind_sl;ind_sl+number_nodes] #position of slack in J
        del_PV = collect(ind_PV.+number_nodes) #position of PV nodes in J

        #delete those positions
        J = J[setdiff(1:end, [del_sl;del_PV]), setdiff(1:end, [del_sl;del_PV])];

        res = J \ [ΔP; ΔQ];

        #With Iwamoto mulipliers the load flow is more robust and does not diverge
        if iwamoto
            iwa = CalcIwamotoMultiplier(J,res,ΔP,ΔQ,Ykk,ind_sl,ind_PV,ind_PQ,number_nodes);
            res *= iwa
        end

        #split result in angle and magnitude
        Δδ  = res[1:length(ΔP)]
        ΔU  = res[length(ΔP)+1:end]

        #Update angle and voltage magnitude, where the latter is the square
        δ[1:end .!= ind_sl, :] += Δδ
        U[setdiff(1:end, [ind_sl; ind_PV]), :] = U[setdiff(1:end, [ind_sl; ind_PV]), :] + ΔU.*U[setdiff(1:end, [ind_sl; ind_PV]), :]

        #Check if reactive power limit is reached; change PV to PQ node
        if Qmax != -1 && Qmin != -1 #Qmax and Qmin are normally vectors with the limits
            #In future, a routine is needed which will get the information about
            #the reactive power limits from each node directly
            if mod(iter,1) == 0 && !isempty(ind_PV_or) #dont change every iteration and check if PV nodes exist
                for i in ind_PV_or # Calc current reactive power at each node
                    Qn[i] = sum(U[i].*U.*Ykk_abs[:,i].*sin.(δ[i].-δ.-θ[:,i]))
                    if Qn[i] >= Qmax[i]
                        index = findall(x-> x==i ,ind_PV) #find the PV node index
                        ind_PV = ind_PV[1:end .!= index, 1]; #delete it from the PV node index list
                        append!(ind_PQ,i)  #append it to the PQ list
                        S_node_gen[i] = real(S_node_gen[i]) + 1im*Qmax[i]
                    elseif Qn[i] <= Qmin[i] #same here of lower limits
                        index = findall(x-> x==i ,ind_PV)
                        ind_PV = ind_PV[1:end .!= index, 1];
                        append!(ind_PQ,i)
                        S_node_gen[i] = real(S_node_gen[i]) + 1im*Qmin[i]
                    elseif isempty(findall(x-> x==i ,ind_PV)) #if inside limits and not in list, set as PV node again
                        index = findall(x-> x==i ,ind_PQ)
                        ind_PQ = ind_PQ[1:end .!= index, 1]; #delete PQ node
                        append!(ind_PV,i)
                    end
                end
            end
        end
    end
end

function NodalAdmittanceMatrice(pg::PowerGrid,U_r_nodes,Ubase)
    Fourpoles = PiModel.(values(pg.lines));
    #changing sign convention from PowerDynamics
    for i in 1:length(Fourpoles)
        Fourpoles[i][1,:] *= -1
    end
    B = Array{Complex{Float64},2}(undef,0,0);
    for i in Fourpoles
        B = vcat(hcat(B,zeros(size(B)[1],2)),hcat(zeros(2,size(B)[1]),i))
    end
    inci = incidence_matrix(powergrid.graph, oriented = false);
    inci_new = zeros(size(inci)[1],2*size(inci)[2]);
    for i in 1:size(inci)[2]
        ind = findall(x->x==1,inci[:,i])
        inci_new[ind[1],2*i-1] = 1
        inci_new[ind[2],2*i] = 1
    end
    #Create voltage ratio matrice Yü that relates the impedances of the different voltage levels
    Yü = zeros(length(U_r_nodes),length(U_r_nodes))
    for i in 1:length(U_r_nodes)
        Yü[:,i] = U_r_nodes[i].*U_r_nodes./(Ubase^2)
    end
    Ykk = Yü.*(inci_new*B*inci_new')
end

function CalcδStartValues(pg,Ykk,ind_sl)
    # Calculates voltage angles based on an 'idle grid' (leerlaufendes netz)
    # Increases power flow robustness for grids with high line impedances
    # more infos (on p. 338):
    # Bernd R. Oswald,
    # Berechnung von Drehstromnetzen, 3. Auflage, Springer Verlag, 2017
    u_sl = 1.
    try #if slack is SlackAlgebraic
        u_sl = collect(values(pg.nodes))[ind_sl].U
    catch #otherwise take slack voltage as 1, could be improved in the future
        u_sl = 1.
    end
    Ykkr  = copy(Ykk)
    Yk_sl = Ykk[:,ind_sl]
    Ykkr[:,ind_sl] .= 0.
    Ykkr[ind_sl,:] .= 0.
    Ykkr[ind_sl,ind_sl] = -Ykk[ind_sl,ind_sl]
    uk   = -inv(Ykkr)*Yk_sl*u_sl
    return angle.(uk)
end

function CalculatePolarLoadFlowJacobian(U,δ,Ykk)
    number_nodes = length(U)
    Ykk_abs = abs.(Ykk)
    θ = angle.(Ykk)
    dPdδ = zeros(number_nodes,number_nodes)
    dPdU = zeros(number_nodes,number_nodes)
    dQdδ = zeros(number_nodes,number_nodes)
    dQdU = zeros(number_nodes,number_nodes)
    for i in 1:number_nodes
        ind_= [collect(1:i-1); collect(i+1:number_nodes)]; #index list without current index

        dPdδ[i,:]  = U[i].*U.*Ykk_abs[:,i].*sin.(δ[i].-δ.-θ[:,i]);
        dPdδ[i,i]  = -1*sum(dPdδ[i,ind_]);

        dPdU[i,:]  = -U[i].*U.*Ykk_abs[:,i].*cos.(δ[i].-δ.-θ[:,i]);
        dPdU[i,i]  = -sum(dPdU[i,ind_]);

        dQdδ[i,:]  = U[i].*Ykk_abs[:,i].*cos.(δ[i].-δ.-θ[:,i]);
        dQdδ[i,i]  = 2.0*U[i].*Ykk_abs[i,i].*cos.(θ[i,i]).+sum(U[ind_].*Ykk_abs[ind_,i].*cos.(δ[i].-δ[ind_].-θ[ind_,i]));

        dQdU[i,:]  = U[i].*Ykk_abs[:,i].*sin.(δ[i].-δ.-θ[:,i]);
        dQdU[i,i]  = (2.0*(U[i].^2).*Ykk_abs[i,i].*sin(-θ[i,i]).-dPdδ[i,i])./U[i];
    end
    return J = [dPdδ dQdδ; dPdU dQdU]
end

function CalcIwamotoMultiplier(J,res,ΔP,ΔQ,Ykk,ind_sl,ind_PV,ind_PQ,number_nodes)
    # For more information:
    #  S. Iwamoto and Y. Tamura,
    # "A Load Flow Calculation Method for Ill-Conditioned Power Systems,"
    # in IEEE Transactions on Power Apparatus and Systems, vol. PAS-100, no. 4,
    # pp. 1736-1743, April 1981, doi: 10.1109/TPAS.1981.316511.
    ΔU = zeros(number_nodes,1)
    ΔU[setdiff(1:end, [ind_sl; ind_PV]), :] = res[length(ΔP)+1:end]

    Δδ = zeros(number_nodes,1)
    Δδ[setdiff(1:end, [ind_sl]), :] = res[1:length(ΔP)]

    ΔUc = ΔU.*exp.(1im*Δδ)
    ΔS  = ΔUc.*(conj.(Ykk)*conj.(ΔUc))
    ΔS  = ΔS[1:end .!= ind_sl, :];

    ΔQ_tmp = zeros(number_nodes,1)
    ΔQ_tmp[intersect(1:end, ind_PQ), :] = ΔQ
    ΔQ_tmp  = ΔQ_tmp[1:end .!= ind_sl, :];

    a = [ΔP;ΔQ_tmp]
    b = -a;
    c = -[real(ΔS);imag(ΔS)]

    g0 = sum(a.*b);
    g1 = sum(b.*b+ 2.0*a.*c);
    g2 = sum(3.0*b.*c);
    g3 = sum(2.0*c.*c);

    f(x) = g3*x^3+g2*x^3+g1*x+g0
    return  find_zero(f,1.0)
end
