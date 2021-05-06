using DifferentialEquations
using ModelingToolkit
using Plots

@variables t, x(t),y(t)
@variables f(x,y)

eqs = f ~ -(cos(x)^2+cos(y)^2)^2

ode = ODESystem(eqs,t,[x,y],[])

rhs = Num(equations(ode)[1].rhs)
grad = Array{Num}(undef,2,1)
ab = Differential.(states(ode))
for (in,dx) in enumerate(ab)
    grad[in] = Num.(expand_derivatives.(map(dx,rhs)))
end

function Substitute(syms::Array{Num},subs_args::Array{Pair{Num,Float64},1}) #SymbolicUtils.Symbolic{Real}
   return Symbolics.value.(substitute.(syms,(subs_args,)))
end

x_range = -80.0:5.:80
y_range = -80.0:5.:80

function test(x_range,y_range,grad)
   tmp = plot()
   for (ind1,i) in enumerate(y_range)
      for (ind2,j) in enumerate(x_range)
         xy = Substitute(grad,[x=>i;y=>j])
         #phase_portrait[ind1,ind2] = (x1[1],x1[2])
         x1 = Float64(xy[1])
         y1 = Float64(xy[2])
         V  = abs(y1/x1)
         dx = sign(x1)*1.25/(sqrt(1.0+V^2))
         dy = 0.0;
         isinf(V) ? dy = 1.25 : dy = V*dx*sign(y1)
         #@assert x1 != 0 "und dy1 = $(dy); dx = $(dx); i = $i, j = $j"
         plot!(f,[j-dx/2.0,j+dx/2.0],[i-dy/2.0,i+dy/2.0]) #,arrow =true
      end
   end
   tmp
end
z = test(x_range,y_range,grad);
z.o
#quiver!(x1-dx/2.0,y-dy/2.0,quiver=(dx,dy))

gr()
N = 2
xa = rand(1:10,N)
ya = rand(1:10,N)
u = rand(N)
v = rand(N)
scatter(xa,ya)
quiver!(xa,ya,quiver=(u,v),linewidth = sqrt.(u.^2+v.^2))
plot!(f,[1.1,2.3],[1.1,2.3], arrow=true)

function a()
   f = plot()
   for i in 1:5
      for j in 1:4
         plot!(f,[i,i+1],[j,j-2],arrow =true,legend=false)
      end
   end
   f
end
a()
