using NLsolve

function sat(x1)
    S10  = 0.1
    S12 = 0.3
    A = (1.2 - sqrt(1.2*S12/S10)) / (1.0 - sqrt(1.2*S12/S10))
    B = S10 / (1 - A)^2
    if x1 > A
        return B*(x1-A)^2/x1
    else
        return 0.0
    end
end

function CalcFieldVoltage(Vrterm,Viterm,Irterm,Iiterm)
    function f!(F,x)
        Ra = 0.0
        Xd = 2.2
        Xq = 2.0
        Xdss = 0.2
        Xqss = 0.2
        Xl = 0.15
        Sd = sat(sqrt(x[1]^2 + x[2]^2))
        Sq = (Xq)/(Xd) * Sd
        Xdsatss = (Xdss - Xl) /  (1 + Sd) + Xl
        Xqsatss = (Xqss - Xl) /  (1 + Sq) + Xl

        #Formel soll null ergeben
        F[1] = Vrterm + Ra*Ir - Xqsatss * Ii - x[1] #x[1] entspricht Vr oder Ed''
        F[2] = Viterm + Ra*Ii + Xdsatss * Ir - x[2] #x[1] entspricht Vi oder Eq''
    end
    sol = nlsolve(f!, [0.0,0.0],ftol = 1e-12)
    return sol.zero
end

Vrterm = 0.99962126337 # aus Lastfluss
Viterm = -0.027519627414 # aus Lastfluss
Ir = 1.0101502178 # aus Lastfluss
Ii = -0.60493032598 # aus Lastfluss

Vr,Vi = CalcFieldVoltage(Vrterm,Viterm,Ir,Ir)

#Berechnung des Polradwinkels
Xq = 2.0
Xd = 2.2
Xqss = 0.2
K = (Xq-Xqss)/(1 + sat(sqrt(Vr^2 + Vi^2)))

δ = atan((Vi+Ir*K)/(Vr-Ii*K))*180/pi

\
using CSV
using DataFrames
using Plots
data = Matrix(CSV.read("C:\\Users\\liemann\\Desktop\\Test_EMT.csv",DataFrame))

t = data[:,1]
uabc = data[:,5:7]
δ = 100*pi.*t


function dq0(a,b,c,δ)
    dqz = sqrt(2/3).*[cos.(δ) cos.(δ.-2/3*pi) cos.(δ.+2/3*pi);
             -sin.(δ) -sin.(δ.-2/3*pi) -sin.(δ.+2/3*pi);
             0.5 0.5 0.5 ] * [a;b;c]
    return dqz
end

tmp = dq0.(uabc[:,1],uabc[:,2],uabc[:,3],δ)
tmp2 = zeros(length(tmp),3)
for (ind,val) in enumerate(tmp)
    tmp2[ind,1] = val[1]
    tmp2[ind,2] = val[2]
    tmp2[ind,3] = val[3]
end


plot(t,sqrt.(tmp2[:,1].^2 +tmp2[:,2].^2))
xlims!((0.099,0.15))

tmp2 = [(1,2); (3,4)]
