

function GetCallbacks(pg::PowerGrid)
    cb = CallbackSet()
    for (ind,val) in enumerate(pg.nodes)
        if typeof(val[2]) == ThreePhaseFault
            cb_tmp = Callback(val[2])
            cb = CallbackSet(cb,cb_tmp)
        end
    end
    return cb
end

function Callback(node::ThreePhaseFault)

    t_on = 1.0
    t_off = 2.0
    
    fun1 = "function timer_hit(u,t,integrator)
                t == $t_on
            end"

    fun2 = "function timer_hit(u,t,integrator)
                t == $t_off
            end"

    fun3 = "function Shortcircuiton(integrator)
                display(1.0)
                #integrator.p[$(node.p_ind[1])] = 80.0
                #integrator.p[$(node.p_ind[2])] = 80.0
                #initialize_dae!(integrator,BrownFullBasicInit())
            end"

    fun4 = "function Shortcircuitoff(integrator)
                display(2.0)
                #integrator.p[$(node.p_ind[1])] = 8e3
                #integrator.p[$(node.p_ind[2])] = 8e3
                #initialize_dae!(integrator,BrownFullBasicInit())
            end"

    cb1 = ContinuousCallback(eval(Meta.parse(fun1)),eval(Meta.parse(fun3)))
    cb2 = ContinuousCallback(eval(Meta.parse(fun2)),eval(Meta.parse(fun4)))
    return CallbackSet(cb1,cb2)
end
