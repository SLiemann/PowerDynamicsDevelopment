module TestModule
    export f1
    export f2

    function f1(a,b)
        c = f2(b)
        return a+b+c
    end

    function f2(a)
        return a*a
    end
end
