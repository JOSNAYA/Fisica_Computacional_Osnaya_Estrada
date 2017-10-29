
__precompile__() # Este comando es para que julia precompile el paquete
module milibreria

using SymPy

export NewtonSympy

"""Método de Newton, se tiene 3 argumentos el primero es la función,xo es una condicion inicial o un arreglo de condiciones iniciales, y tolerancia para el error permitido del método."""

function Newton(f,xo,tolerancia) #Función, valor inicial, intervalo
ERROR=[]
x=Sym("x")
df=lambdify(diff(f(x),x))
   while abs(f(xo))>tolerancia
        xo=xo-f(xo)/df(xo)
        push!(Xk,xo)
    end    
    xo
end

export Integraltrapecio
 """Método del trapecio, a-punto inicial a-punto final del intervalo,n nos dice en cuantas partes se va a dividir el intervalo.""" 
function Integral_trapecio(f,a,b,n)
h=(b-a)/n
valor=h*(f(a)+f(b))/2   
    for i in 1:n
      valor=valor+f(a+((b-a)*i/n))*h
    end
    valor
end

export Integralsimpson
"""Método de Simpson de segundo orden,el primero es la función, a-punto inicial b-punto final del intervalo, n-cuantas partes se va a dividir el intervalo. """
function Integral_simpson(f,a,b,n)
h=(b-a)/n
intervalo=linspace(a,b,n+1)    
s=convert(Int64,n/2)  
valor=0    
    for i in 1:s
        valor=valor+(f(intervalo[2i-1])+4f(intervalo[2i])+f(intervalo[2i+1]))*(h/3)
    end
    valor
end 

export Integralriemann
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

export derivada_orden

""" Derivada numérica con error de orden h^n. El método admite una funcion f(x), un punto x donde se calcula la derivada, una diferencia h con la que se aproxima la derivada, y orden es el orden de error del método como potencia de h."""
function coeficiente_taylor_asinh(i)
x=Sym("x");
h=Sym("h");
N(subs(diff(asinh(x/2),x,i)/factorial(Float64(i)),x=>0))
end

function diferencia_finita_simetrica(f)
    return (x,h)->f(x+h,h)-f(x-h,h)
end

function derivada_orden(f,x,h,orden)
    F(x,h)=f(x)
    df=0
    for i in 1:orden
        F=diferencia_finita_simetrica(F)
        df=df+coeficiente_taylor_asinh(i)*F(x,h)
    end
    return df/h
end

export MetodoEuler

""" Método de Euler idependiente de la dimensión para resolver una ecuación diferencial. El primero es la función f, condiciones inicial, por último h es el paso del método."""

function MetodoEuler(f,xo,to,t,h)
    listt=to:h:t
    listx=[xo]
    x=xo
    for i in 1:length(listt)-1
    x=x+h*f(x,listt[i])
        push!(listx,x)
    end
    return listt,listx
end

export MetodoEuler_implicito
""" Método de Euler idependiente para resolver una ecuación diferencial en un intervalo (to,t). El método admite una función f, y una condicion inicial xo, por último h es el paso del método. 
Como resultado se tiene un arreglo con dos entradas, la primera entrada es el arreglo de valores de t, la segunda entra es un arreglo con los valores de x"""

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


export RungeKutta4th
""" Método de Runge Kutta a cuarto orden. El método puede recibir un arreglo de funciones f, un arreglo de condiciones iniciales xo, un punto inicial y final to y t correspondientemente y un paso h.
Como resultado se tiene un arreglo con dos entradas, la primera entrada es el arreglo de valores de t, la segunda entra es un arreglo con arreglos cada uno con el valor de cada paso dado para cada variable"""

function RungeKutta4th(f,xo,to,t,h)
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


