module Quadrature

function quad_hammer(f,x₀::Vector,y₀::Vector,degree::Int=2)
    if degree==1
        a=[1/3,1/3,1/3]
        x=x₀'*a
        y=y₀'*a
        return f(x,y)/2
    elseif degree==2
        a=[2/3,1/6,1/6]
        b=[1/6,2/3,1/6]
        c=[1/6,1/6,2/3]
        w=[1/3,1/3,1/3]
        x=x₀'*[a b c]
        y=y₀'*[a b c]
        return reduce(+,w'.*f.(x,y))/2
    elseif degree==3
        a=[1/3,1/3,1/3]
        b=[0.6,0.2,0.2]
        c=[0.2,0.6,0.2]
        d=[0.2,0.2,0.6]
        w=[-27/48,25/48,25/48,25/48]
        x=x₀'*[a b c d]
        y=y₀'*[a b c d]
        return reduce(+,w'.*f.(x,y))/2
    end
end

end
