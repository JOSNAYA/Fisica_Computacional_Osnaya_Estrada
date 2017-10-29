
__precompile__() 

module herramientas

using SymPy

export Newton
"""Método de Newton,función,xo es una condicion inicial, y tolerancia para el error permitido del método."""

function Newton(f,xo,tolerancia) 
valor=[]
x=Sym("x")
df=lambdify(diff(f(x),x))
   while abs(f(xo))>tolerancia
        xo=xo-f(xo)/df(xo)
        push!(valor,xo)
    end    
    xo
end

export Integral_trapecio
 """Método del trapecio, a-punto inicial a-punto final del intervalo,n nos dice en cuantas partes se va a dividir el intervalo.""" 
function Integral_trapecio(f,a,b,n)
h=(b-a)/n
valor=h*(f(a)+f(b))/2   
    for i in 1:n
      valor=valor+f(a+((b-a)*i/n))*h
    end
    valor
end

export Integral_simpson
"""Método de Simpson de segundo orden,el primero es la función, a-punto inicial b-punto final del intervalo, n-cuantas partes se va a dividir el intervalo. """
function Integral_simpson(f,a,b,n)
valor=0.0 
valor=((b-a)/6) * (f(a)+4(f((a+b)/2))+f(b))
return valor
end 

export Integral_riemann
"""Método de integración de Riemann, el primero es la función, a-punto inicial b-punto final del intervalo, n-partes se va a dividir el intervalo. """
function Integral_riemann(f,a,b,n)
intervalo=linspace(a,b,n)
valor=0
s=convert(Int64,n-1)    
    for i in 1:s
        valor=valor+f(intervalo[i+1])*(intervalo[i+1]-intervalo[i])
    end
    valor
end

export Euler

""" Método de Euler idependiente de la dimensión. El primero es la función f, condiciones inicial, h es el paso del método."""

function Euler(f,xo,to,t,h)
    listt=to:h:t
    listx=[xo]
    x=xo
    for i in 1:length(listt)-1
    x=x+h*f(x,listt[i])
        push!(listx,x)
    end
    return listt,listx
end

export Euler_implicito
""" Método de Euler idependiente de la dimensión. El primero es la funcion, condicion inicial xo, h es el paso del método.""" 

function Newton_para_implicito(f,x0,t) 
listx=[]
df(y,t)=N(diff(f(x,t),x))              
    while abs(f(x0,t))>0.000001
        x0=x0-f(x0,t)/df(x0,t)         
        push!(listx,x0)                
    end    
    x0
end

function Euler_implicito(f,x0,t0,t,h)
    listt=t0:h:t                         
    listx=zeros(length(listt))           
    listx[1] = x0                        
    for i in 2:length(listt)
        g(x,s)=x-h*f(x,listt[i-1])-x0    
        x0=x0+h*f(Newton(g,x0,listt[i-1]),listt[i])     
        listx[i]=x0
    end
    return listt,listx
end


export Runge_Kutta_4
""" Método de Runge Kutta de cuarto orden. El primero es funcion , condicion inicial, un punto inicial y un punto final y  h."""

function Runge_Kutta_4(f,xo,to,t,h)
    x=xo
    listt=to:h:t
    listx=[xo]
    for i in 1:length(listt)-1
        
    k1=f(x,listt[i])
    k2=f(x+(h/2)*k1,listt[i]+h/2) 
    k3=f(x+(h/2)*k2,listt[i]+h/2) 
    k4=f(x+h*k3,listt[i]+h)
        
    x=x+(h/6)*(k1+2k2+2k3+k4)
        
        push!(listx,x)    
    end
    return listt,listx
end

end
