using PowerDynamics
using PowerDynamics: guess
using LightGraphs #incidence_matrix
using Roots #for iwamoto_multiplier
using LinearAlgebra

#Pi models for nodal admittance matrice
PiModel(L::PiModelLine) = PiModel(L.y,L.y_shunt_km,L.y_shunt_mk,1,1)
PiModel(T::Transformer) = PiModel(T.y,0,0,T.t_ratio,1)
PiModel(S::StaticLine)  = PiModel(S.Y,0,0,1,1)
PiModel(R::RLLine)      = PiModel(1/(R.R+R.ω0*R.L),0,0,1,1)

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

PowerNodeLoad(S::SlackAlgebraic,U) = 0.  #treated as generation
PowerNodeLoad(F::FourthOrderEq,U) = 0.
PowerNodeLoad(F::FourthOrderEqExciterIEEEDC1A,U)  = 0.
PowerNodeLoad(F::FourthOrderEqGovernorExciterAVR,U)  = 0.
PowerNodeLoad(F::FourthOrderEqGovernorIEEEG1,U)  = 0.
PowerNodeLoad(S::SwingEq,U)  = 0.
PowerNodeLoad(V::VSIMinimal,U) = complex(V.P,(abs(U)-V.V_r)/V.K_Q+V.Q)   #treated as load with changed sign to include voltage dependency
PowerNodeLoad(V::VSIVoltagePT1,U) = complex(V.P,(abs(U)-V.V_r)/V.K_Q+V.Q)   #treated as load with changed sign to include voltage dependency
PowerNodeLoad(L::PVAlgebraic,U) = -L.P  #treated as load
PowerNodeLoad(L::PQAlgebraic,U) = -complex(L.P,L.Q) #treated as load
PowerNodeLoad(L::VoltageDependentLoad,U) = -complex(L.P, L.Q) * (L.A * abs(U)^2 + L.B * abs(U) + 1 - L.A - L.B)
PowerNodeLoad(L::ExponentialRecoveryLoad,U)  = -(L.P0*((abs(U)/L.V0)^L.Nps) + 1im*L.Q0*((abs(U)/L.V0)^L.Nqs))
PowerNodeLoad(L::CSIMinimal,U)  = -U*conj(L.I_r)

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
PowerNodeGeneration(L::VoltageDependentLoad) = 0.
PowerNodeGeneration(L::ExponentialRecoveryLoad)  = 0.
PowerNodeGeneration(L::CSIMinimal,U)  = 0.


function PowerFlowClassic(pg::PowerGrid, U_r_nodes::Vector{Float64}; Ubase=380e3, ind_sl::Int64 = 0,max_tol = 1e-6,iter_max::Int64  = 30, Qmax=-1., iwamoto::Bool =false)
    nodetypes = NodeType.(values(pg.nodes))
    if !isempty(findall(x-> x==0, nodetypes))
        ind_sl = findall(x-> x==0, nodetypes)[1]
        @info "Reference node is bus no. $ind_sl"
    end
    ind_PV_or = findall(x-> x==1, nodetypes)
    ind_PQ_or = findall(x-> x==2, nodetypes)

    ind_PQ = ind_PQ_or #copy is needed for checking Qlimit
    ind_PV = ind_PV_or #copy is needed for checking Qlimit

    PQkeys = findall( x-> typeof(x) == VoltageDependentLoad, pg.nodes)
    PVkeys = findall( x-> typeof(x) == FourthOrderEq , pg.nodes)
    append!(PVkeys, findall( x-> typeof(x) == FourthOrderEqreal , pg.nodes))

    number_nodes = length(pg.nodes);
    U = ones(number_nodes,1);
    δ = zeros(number_nodes,1);
    Ykk = NodalAdmittanceMatrice(pg,U_r_nodes,Ubase);
    Ykk_abs = abs.(Ykk);
    θ   = angle.(Ykk);
    Pn  = similar(U);
    Qn  = similar(U);
    S_node_gen = Array{Complex{Float64},2}(undef,number_nodes,1)
    S_node_load = Array{Complex{Float64},2}(undef,number_nodes,1)
    #is called only once, since it should only changes, when power limits are reached
    for (ind,val) in enumerate(values(pg.nodes)) S_node_gen[ind] = PowerNodeGeneration(val) end

    #Start der Iterationen
    for iter in 1:iter_max+1
        for i in 1:number_nodes # Knotenleistungen berechnen, anhand der aktuellen Spannung
            Pn[i] = sum(U[i].*U.*Ykk_abs[:,i].*cos.(δ[i].-δ.-θ[:,i]))
            Qn[i] = sum(U[i].*U.*Ykk_abs[:,i].*sin.(δ[i].-δ.-θ[:,i]))
        end
        #get current load powers
        for (ind,val) in enumerate(values(pg.nodes)) S_node_load[ind] = PowerNodeLoad(val,U[ind]*exp(1im*δ[ind])) end
        S_node = S_node_gen + S_node_load

        ΔP = real(S_node) - Pn
        ΔQ = imag(S_node) - Qn

        ΔP = ΔP[1:end .!= ind_sl, :];
        ΔQ = ΔQ[setdiff(1:end, [ind_sl; ind_PV]), :];

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
        del_sl = [ind_sl;ind_sl+number_nodes] #Position des Slacks in J
        del_PV = collect(ind_PV.+number_nodes) #Position der PV Knoten in J

        J = J[setdiff(1:end, [del_sl;del_PV]), setdiff(1:end, [del_sl;del_PV])];

        res = J \ [ΔP; ΔQ];

        iwamoto_multiplier = ones(length([ΔP; ΔQ]),1)
        if iwamoto
            iwa = CalcIwamotoMultiplier(J,res,ΔP,ΔQ,Ykk,ind_sl,ind_PV,ind_PQ,number_nodes);
            res *= iwa
        end

        Δδ  = res[1:length(ΔP)] #more elegant way, like (ΔU,Δδ) .= J \ [ΔP; ΔQ] ?
        ΔU  = res[length(ΔP)+1:end]

        δ[1:end .!= ind_sl, :] += Δδ
        U[setdiff(1:end, [ind_sl; ind_PV]), :] = U[setdiff(1:end, [ind_sl; ind_PV]), :] + ΔU.*U[setdiff(1:end, [ind_sl; ind_PV]), :]

        #Funktion raus machen
        #Check if reactiv limit is achieved; change PV to PQ node
        Qmax =  ones(number_nodes,1)*Inf #dummy
        Qmin = -ones(number_nodes,1)*Inf #dummy
        if mod(iter,2) == 0 && !isempty(ind_PV_or) #dont change every iter and check if PV nodes exist
            for i in ind_PV_or # Knotenleistungen berechnen, anhand der aktuellen Spannung
                Qn[i] = sum(U[i].*U.*Ykk_abs[:,i].*sin.(δ[i].-δ.-θ[:,i]))
                # check both limits of Qmax : [-Inf ]
                if Qn[i] >= Qmax[i]
                    index = findall(x-> x==i ,ind_PV)
                    ind_PV = ind_PV[1:end .!= index, 1];
                    append!(ind_PQ,i)
                    S_node_gen[i] = real(S_node_gen[i]) + 1im*Qmax[i] #VIELLEICHT IST DIESE ART DES LÖSCHENS FALSCH!!!!!!
                elseif Qn[i] <= Qmin[i]
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

function NodalAdmittanceMatrice(pg::PowerGrid,U_r_nodes,Ubase)
    Fourpoles = PiModel.(values(pg.lines)); # get all fourpoles
    #changing sign convention
    for i in 1:length(Fourpoles)
        Fourpoles[i][1,:] *= -1
    end
    B = Array{Complex{Float64},2}(undef,0,0);
    for i in Fourpoles  #create diagonal matrice of fourpoles
        B = vcat(hcat(B,zeros(size(B)[1],2)),hcat(zeros(2,size(B)[1]),i))
    end
    inci = incidence_matrix(powergrid.graph, oriented = false);
    inci_new = zeros(size(inci)[1],2*size(inci)[2]);
    for i in 1:size(inci)[2]  #change incidence matrice so that the sum of each column is one
        ind = findall(x->x==1,inci[:,i])
        inci_new[ind[1],2*i-1] = 1
        inci_new[ind[2],2*i] = 1
    end
    #Create matrice that transforms the different voltage levels
    Yü = zeros(length(U_r_nodes),length(U_r_nodes))
    for i in 1:length(U_r_nodes)
        Yü[:,i] = U_r_nodes[i].*U_r_nodes./(Ubase^2)
    end
    Ykk = Yü.*(inci_new*B*inci_new')
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
        ind_= [collect(1:i-1); collect(i+1:number_nodes)];

        dPdδ[i,:]  = U[i].*U.*Ykk_abs[:,i].*sin.(δ[i].-δ.-θ[:,i]);
        dPdδ[i,i] = -1*sum(dPdδ[i,ind_]);

        dPdU[i,:]  = -U[i].*U.*Ykk_abs[:,i].*cos.(δ[i].-δ.-θ[:,i]);
        dPdU[i,i] = -sum(dPdU[i,ind_]);

        dQdδ[i,:]  = U[i].*Ykk_abs[:,i].*cos.(δ[i].-δ.-θ[:,i]);
        dQdδ[i,i] = 2.0*U[i].*Ykk_abs[i,i].*cos.(θ[i,i]).+sum(U[ind_].*Ykk_abs[ind_,i].*cos.(δ[i].-δ[ind_].-θ[ind_,i]));

        dQdU[i,:]  = U[i].*Ykk_abs[:,i].*sin.(δ[i].-δ.-θ[:,i]);
        dQdU[i,i] = (2.0*(U[i].^2).*Ykk_abs[i,i].*sin(-θ[i,i]).-dPdδ[i,i])./U[i];
    end
    return J = [dPdδ dQdδ; dPdU dQdU]
end

function CalcIwamotoMultiplier(J,res,ΔP,ΔQ,Ykk,ind_sl,ind_PV,ind_PQ,number_nodes)
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
