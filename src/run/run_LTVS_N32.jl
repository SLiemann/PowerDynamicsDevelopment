using PowerDynamics
using Plots
using IfElse

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
end
begin
    pg = LTVS_Test_System_N32()
    Qmax   = [Inf, Inf, Inf,Inf,Inf*53.0*sqrt(1-0.8377358^2),Inf]
    Qmin   = -Qmax
    U,δ,ic0 = PowerFlowClassic(pg,iwamoto = false,max_tol = 1e-4,iter_max = 30,Qmax = Qmax, Qmin = Qmin,Qlimit_iter_check=3)
    #pg, ic0 = InitializeInternalDynamics(pg,ic0)
    #pgsol  = run_LTVS_N32_simulation(pg,ic0,(0.0,300.0))
end

CalcEigenValues(pg,[],plot =true, output = true)
xlims!((-3,0))
plot(pgsol,["bus2","bus3","bus4","bus5"],:v, legend = legend=:bottomright)
plot(pgsol,"bus5",:timer)
plot(pgsol,"bus5",:ifd)


Uc = U.*exp.(1im*δ/180*pi)
Ykk = NodalAdmittanceMatrice(pg)
S  = round.(Uc.*(conj.(Ykk)*conj.(Uc)),digits=5)
abs.(S)


cd(joinpath(DEPOT_PATH[1], "registries", "General.toml")) do
           deps = Pkg.dependencies()
           registry = Pkg.TOML.parse(read("Registry.toml", String))
           general_pkgs = registry["packages"]

           constrained = Dict{String, Tuple{VersionNumber,VersionNumber}}()
           for (uuid, dep) in deps
               suuid = string(uuid)
               dep.is_direct_dep || continue
               dep.version === nothing && continue
               haskey(general_pkgs, suuid) || continue
               pkg_meta = general_pkgs[suuid]
               pkg_path = joinpath(pkg_meta["path"], "Versions.toml")
               versions = Pkg.TOML.parse(read(pkg_path, String))
               newest = maximum(VersionNumber.(keys(versions)))
               if newest > dep.version
                   constrained[dep.name] = (dep.version, newest)
               end
           end

           return constrained
end
