using NonlinearSolve, Plots

function f_Zeq(u,p)
    Rload = u[1]
    Xload = p[1]
    Rl = p[2]
    Xl = p[3]
    Xcf = p[4]
    Xlf = p[5]
    Rf = p[6]

    Zgrid = (Rload+1im*Xload)/(Rload*1im*Xload) + Rl + 1im*Xlf
    Zeq = (Zgrid + -1im*Xcf)/(-1im*Xcf*Zgrid) + Rf +1im*Xlf
    1.0 - abs(Zeq)
end


RX = 0.1
Zline = 0.05
Xlv = Zline/sqrt(1+RX^2)
Rlv = RX*Xlv

Xcfv = 5.3052
Xlfv = 0.0314
Rfv = 5.0000e-04
p = [1e90 Rlv Xlv Xcfv Xlfv Rfv]
u0 = 1.0
probN = NonlinearProblem(f_Zeq, u0, p,reltol=1e-9)
sol = solve(probN, NewtonRaphson())

Rloadmin = 0.987114216244142 # wenn Xload nicht vorhanden
Xloadmin = 0.76145901778670 # wenn Rload nicht vorhanden 

r = []
x = []
for Xloadv =4:-0.0001:0.78
    u0 = 0.95
    p = [Xloadv Rlv Xlv Xcfv Xlfv Rfv]
    probN = NonlinearProblem(f_Zeq, u0, p)
    sol = solve(probN, NewtonRaphson(),reltol=1e-7)
    if ~isnan(sol.u)
     #r = plot!([sol.u],[Xloadv],seriestype=:scatter, legend=false)
        r  = vcat(r,sol.u)
        x  = vcat(x,Xloadv)
    end
end
plot(r,x)



#using ModelingToolkit
#@parameters Rload Xload Rl Xl Xcf Xlf Rf

#Zgrid = (Rload+1im*Xload)/(Rload*1im*Xload) + Rl + 1im*Xlf
#Zeq = (Zgrid + -1im*Xcf)/(-1im*Xcf*Zgrid) + Rf +1im*Xlf
#abs(Zeq)