using PlotlyJS, DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32.jl")
    nothing
end

pg_normal = LTVS_Test_System_N32();
pg_fault = GetPostFaultLTVSPG(pg_normal);
begin
    Qmax   = [Inf,Inf, Inf, Inf,Inf,5300/8000*sqrt(1-0.85^2),Inf]
    Qmin   = -Qmax
    U1,δ1,ic0,cu = PowerFlowClassic(pg_normal,iwamoto = false,max_tol = 1e-4,iter_max = 100,Qmax = Qmax, Qmin = Qmin,Qlimit_iter_check=80)
    pg_dyn, ic0 = InitializeInternalDynamics(pg_normal,ic0)
    pgsol  = run_LTVS_N32_simulation(pg_dyn,ic0,(0.0,165.0));
    or = plotv(pgsol,"bus_load");
end

plot(plotv(pgsol,"bus_load"))
Uv = Vector{Float64}()
bus_ind = 5
Qmax   = [Inf,Inf, Inf, Inf,Inf,5300/8000*sqrt(1-0.85^2),Inf]
Qmin   = -Qmax
for i=6:9
    pg_normal = LTVS_Test_System_N32(tap=i);
    pg_fault = GetPostFaultLTVSPG(pg_normal);
    U1,δ1,ic0,cu = PowerFlowClassic(pg_fault,iwamoto = false,max_tol = 1e-8,iter_max = 100,Qmax = Qmax, Qmin = Qmin,
                                    Qlimit_iter_check=80)
                                    push!(Uv,U1[bus_ind])
end

time = pgsol.dqsol.t;
len = size(time)[1];
pf = Vector{GenericTrace}()
for i in Uv
    tmp = scatter(x=time,y=i.*ones(len))
    push!(pf,tmp)
end
plot([or;pf])





pg_normal = LTVS_Test_System_N32();
pg_fault = GetPostFaultLTVSPG(pg_normal);

Qmax   = [Inf,Inf, Inf, Inf,Inf,5300/8000*sqrt(1-0.85^2),Inf]
Qmin   = -Qmax
U1,δ1,ic0,cu = PowerFlowClassic(pg_fault,iwamoto = false,max_tol = 1e-4,iter_max = 100,Qmax = Qmax, Qmin = Qmin,Qlimit_iter_check=80)
pg_dyn, ic1 = InitializeInternalDynamics(pg_fault,ic0)

ic2 = find_valid_initial_condition(pg_dyn,ic1)

sqrt(ic2[11]^2+ic2[12]^2)
sqrt(ic1[11]^2+ic1[12]^2)


