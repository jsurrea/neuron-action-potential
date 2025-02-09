"""PROYECTO FINAL"""
# Programación Científica
# Universidad de los Andes
# Isai Daniel Chacón Silva   201912015
# Juan Sebastián Urrea López 201914710
# Erick Sebastián Lozano Roa 201913084

# MODELO MATEMÁTICO DEL POTENCIAL DE ACCIÓN
# Librerías

import numpy as np
import struct as st
import scipy.optimize as opt
from scipy.integrate import solve_ivp, odeint

# Variables

paso = 0.01     # Define la longitud del paso de iteración en cada método.

# Valores alfa y beta

def beta_m(V):
    return 4*np.exp(-(V+65)/18)
def alfa_h(V):
    return 0.07*np.exp(-(V+65)/20)
def beta_h(V):
    return 1/(1 + np.exp(-(V+35)/10))
def alfa_n(V):
    return 0.01*(V + 55)/(1 - np.exp(-(V+55)/10))
def beta_n(V):
    return 0.125*np.exp(-(V+65)/80)
def alfa_m(V):
    return 0.1*(V + 40)/(1 - np.exp(-(V+40)/10))

# Ecuaciones diferenciales

def derivative_maker(temperatura, funcion_corriente):
    """
    Devuelve la función de derivadas de las variables v,n,m,h con los parámetros de temperatura y corriente dados.
    :param temperatura: Number. Se usa para calcular el factor de temperatura en los canales iónicos.
    :param funcion_corriente: Function. Debe recibir como parámetro el tiempo y devolver un número.
    :return: Function. Vector de derivadas evaluadas para un tiempo y valor de las variables particular.
    """
    g_Na = 120                              # Conductancia por unidad de área del canal de Na+.
    g_K = 36                                # Conductancia por unidad de área del canal de K+.
    g_L = 0.3                               # Conductancia por unidad de área del canal de fuga.
    E_Na = 50                               # Potencial de Nernst de Na+.
    E_K = -77                               # Potencial de Nernst de K+.
    E_L = -54.4                             # Potencial de Nernst de fuga.
    C_M = 1                                 # Capacitancia de la membrana celular.
    Q10 = 3                                 # Radio de las tasas por un incremento de temperatura de 10ºC.
    Tbase = 6.3                             # Temperatura base.
    phi = Q10**((temperatura - Tbase)/10)   # Factor de temperatura.

    def dVar(t,variables):
        """
        Devuelve un vector con las derivadas de las funciones de V,n,m,h indexadas en ese orden.
        :param t: Number. Tiempo en el que se evaluarán las derivadas.
        :param variables: Iterable. Valor de las variables en el orden v,n,m,h.
        :return: ndarray. Vector de derivadas respecto al tiempo en el orden v,n,m,h.
        """
        v,n,m,h = variables                                                                         # Unpacking.
        vect = np.zeros((4,))                                                                       # Inicializar.
        vect[0] = (funcion_corriente(t)-g_L*(v-E_L)-g_K*(n**4)*(v-E_K)-g_Na*(m**3)*h*(v-E_Na))/C_M  # Vm.
        vect[1] = phi*(alfa_n(v)*(1-n)-beta_n(v)*n)                                                 # N.
        vect[2] = phi*(alfa_m(v)*(1-m)-beta_m(v)*m)                                                 # M.
        vect[3] = phi*(alfa_h(v)*(1-h)-beta_h(v)*h)                                                 # H.
        return vect
    return dVar

# Gráficas

def EulerFor(initial_conditions, tf, temperatura, funcion_corriente):
    """
    Implementa Euler Forward para resolver el sistema de ecuaciones diferenciales.
    :param initial_conditions: Iterable. Condición inicial de las funciones en el orden v,n,m,h.
    :param tf: Number. Tiempo final o rango de tiempo de la simulación.
    :param temperatura: Number. Se usa para calcular el factor de temperatura en los canales iónicos.
    :param funcion_corriente: Function. Debe recibir como parámetro el tiempo y devolver un número.
    :return: Tuple. Primer elemento es el vector de tiempo, segundo elemento es la solución de Vm.
    """
    t = np.arange(0.0, tf+paso, paso)                                           # Crear vector de tiempo.
    dVar = derivative_maker(temperatura, funcion_corriente)                     # Crear función de derivadas.
    variables = np.zeros((4,np.size(t)))                                        # Inicializar.
    variables[:,0] = initial_conditions                                         # Condiciones iniciales.
    for i in range(1, np.size(t)):                                              # Iterar.
        variables[:,i] = variables[:,i-1] + paso*dVar(t[i],variables[:,i-1])    # Ecuación de avance.
    return t, variables[0,:]

def EulerBack(initial_conditions, tf, temperatura, funcion_corriente):
    """
    Implementa Euler Backward para resolver el sistema de ecuaciones diferenciales.
    :param initial_conditions: Iterable. Condición inicial de las funciones en el orden v,n,m,h.
    :param tf: Number. Tiempo final o rango de tiempo de la simulación.
    :param temperatura: Number. Se usa para calcular el factor de temperatura en los canales iónicos.
    :param funcion_corriente: Function. Debe recibir como parámetro el tiempo y devolver un número.
    :return: Tuple. Primer elemento es el vector de tiempo, segundo elemento es la solución de Vm.
    """
    t = np.arange(0.0, tf+paso, paso)                                               # Crear vector de tiempo.
    dVar = derivative_maker(temperatura, funcion_corriente)                         # Crear función de derivadas.
    variables = np.zeros((4,np.size(t)))                                            # Inicializar.
    variables[:,0] = initial_conditions                                             # Condiciones iniciales.
    for i in range(1, np.size(t)):                                                  # Iterar.
        variables[:,i] = opt.fsolve(lambda act,ant,tem:ant+paso*dVar(tem,act)-act,  # Función a resolver con fsolve.
                                    variables[:,i-1], (variables[:,i-1], t[i]))     # Variables y parámetros de fsolve.
    return t, variables[0,:]

def EulerMod(initial_conditions, tf, temperatura, funcion_corriente):
    """
    Implementa Euler Modificado para resolver el sistema de ecuaciones diferenciales.
    :param initial_conditions: Iterable. Condición inicial de las funciones en el orden v,n,m,h.
    :param tf: Number. Tiempo final o rango de tiempo de la simulación.
    :param temperatura: Number. Se usa para calcular el factor de temperatura en los canales iónicos.
    :param funcion_corriente: Function. Debe recibir como parámetro el tiempo y devolver un número.
    :return: Tuple. Primer elemento es el vector de tiempo, segundo elemento es la solución de Vm.
    """
    t = np.arange(0.0, tf+paso, paso)                                               # Crear vector de tiempo.
    dVar = derivative_maker(temperatura, funcion_corriente)                         # Crear función de derivadas.
    variables = np.zeros((4,np.size(t)))                                            # Inicializar.
    variables[:,0] = initial_conditions                                             # Condiciones iniciales.
    for i in range(1, np.size(t)):                                                  # Iterar.
        variables[:,i] = opt.fsolve(lambda act,ant,tem:ant+paso/2*(dVar(tem,act)+dVar(tem,ant))-act,
                                    variables[:,i-1], (variables[:,i-1], t[i]))     # Llamar a fsolve como en back.
    return t, variables[0,:]

def RungeKutta2(initial_conditions, tf, temperatura, funcion_corriente):
    """
    Implementa Runge Kutta de segundo orden para resolver el sistema de ecuaciones diferenciales.
    :param initial_conditions: Iterable. Condición inicial de las funciones en el orden v,n,m,h.
    :param tf: Number. Tiempo final o rango de tiempo de la simulación.
    :param temperatura: Number. Se usa para calcular el factor de temperatura en los canales iónicos.
    :param funcion_corriente: Function. Debe recibir como parámetro el tiempo y devolver un número.
    :return: Tuple. Primer elemento es el vector de tiempo, segundo elemento es la solución de Vm.
    """
    t = np.arange(0.0, tf+paso, paso)                           # Crear vector de tiempo.
    dVar = derivative_maker(temperatura, funcion_corriente)     # Crear función de derivadas.
    variables = np.zeros((4,np.size(t)))                        # Inicializar.
    variables[:,0] = initial_conditions                         # Condiciones iniciales.
    for i in range(1, np.size(t)):                              # Iterar.
        k1 = dVar(t[i],variables[:,i-1])                        # Primer ponderado.
        k2 = dVar(t[i]+paso, variables[:,i-1]+paso*k1)          # Segundo ponderado.
        variables[:,i] = variables[:,i-1] + paso/2*(k1+k2)      # Ecuación de avance.
    return t, variables[0,:]

def RungeKutta4(initial_conditions, tf, temperatura, funcion_corriente):
    """
    Implementa Runge Kutta de cuarto orden para resolver el sistema de ecuaciones diferenciales.
    :param initial_conditions: Iterable. Condición inicial de las funciones en el orden v,n,m,h.
    :param tf: Number. Tiempo final o rango de tiempo de la simulación.
    :param temperatura: Number. Se usa para calcular el factor de temperatura en los canales iónicos.
    :param funcion_corriente: Function. Debe recibir como parámetro el tiempo y devolver un número.
    :return: Tuple. Primer elemento es el vector de tiempo, segundo elemento es la solución de Vm.
    """
    t = np.arange(0.0, tf+paso, paso)                                   # Crear vector de tiempo.
    dVar = derivative_maker(temperatura, funcion_corriente)             # Crear función de derivadas.
    variables = np.zeros((4,np.size(t)))                                # Inicializar.
    variables[:,0] = initial_conditions                                 # Condiciones iniciales.
    for i in range(1, np.size(t)):                                      # Iterar.
        k1 = dVar(t[i],variables[:,i-1])                                # Primer ponderado.
        k2 = dVar(t[i]+paso/2, variables[:,i-1]+paso/2*k1)              # Segundo ponderado.
        k3 = dVar(t[i]+paso/2, variables[:,i-1]+paso/2*k2)              # Tercer ponderado.
        k4 = dVar(t[i]+paso, variables[:,i-1]+paso*k3)                  # Cuarto ponderado.
        variables[:,i] = variables[:,i-1] + paso/6*(k1+2*k2+2*k3+k4)    # Ecuación de avance.
    return t, variables[0,:]

def scipy_solve(initial_conditions, tf, temperatura, funcion_corriente):
    """
    Implementa la función solve_ivp de Scipy para resolver el sistema de ecuaciones diferenciales.
    :param initial_conditions: Iterable. Condición inicial de las funciones en el orden v,n,m,h.
    :param tf: Number. Tiempo final o rango de tiempo de la simulación.
    :param temperatura: Number. Se usa para calcular el factor de temperatura en los canales iónicos.
    :param funcion_corriente: Function. Debe recibir como parámetro el tiempo y devolver un número.
    :return: Tuple. Primer elemento es el vector de tiempo, segundo elemento es la solución de Vm.
    """
    dVar = derivative_maker(temperatura, funcion_corriente)                 # Crear función de derivadas.
    sol = solve_ivp(dVar, [0.0,tf], initial_conditions, method='RK45')      # Llamar a solve_ivp de Scipy.
    return sol.t, sol.y[0]

def scipy_odeint(initial_conditions, tf, temperatura, funcion_corriente):
    """
    Implementa la función odeint de Scipy para resolver el sistema de ecuaciones diferenciales.
    :param initial_conditions: Iterable. Condición inicial de las funciones en el orden v,n,m,h.
    :param tf: Number. Tiempo final o rango de tiempo de la simulación.
    :param temperatura: Number. Se usa para calcular el factor de temperatura en los canales iónicos.
    :param funcion_corriente: Function. Debe recibir como parámetro el tiempo y devolver un número.
    :return: Tuple. Primer elemento es el vector de tiempo, segundo elemento es la solución de Vm.
    """
    t = np.arange(0.0, tf+paso, paso)                                       # Crear vector de tiempo.
    dVar = derivative_maker(temperatura, funcion_corriente)                 # Crear función de derivadas.
    sol = odeint(lambda x,tem:dVar(tem,x), initial_conditions, t)[:,0]      # Llamar a odeint de Scipy.
    return t,sol

def exportar(data_figures):
    """Crea archivos binarios de cada gráfica contenida en el diccionario data_figures."""
    for filename,(tiempo,values) in data_figures.items():           # Iterar sobre las funciones y sus datos.
        packing = st.pack("d"*np.size(tiempo),*tiempo)              # Empaquetar vector de tiempo.
        with open("Tiempo_" + filename + ".bin","wb") as file:      # Crear archivo para guardar el tiempo.
            file.write(packing)                                     # Escribir en el archivo.
        packing = st.pack("d"*np.size(values),*values)              # Empaquetar vector de valores.
        with open("Valores_" + filename + ".bin","wb") as file:     # Crear archivo para guardar los valores.
            file.write(packing)                                     # Escribir en el archivo.

def importar(filetiempo, filevalues):
    """
    Importa los datos provenientes de dos archivos binarios codificados en double.
    :param filetiempo: String. Ruta del archivo que contiene los valores del tiempo.
    :param filevalues: String. Ruta del archivo que contiene los valores del potencial.
    :return: Tuple. Primer elemento es el vector de tiempo, segundo elemento es el vector de Vm.
    """
    with open(filetiempo,"rb") as file:                             # Abrir archivo del tiempo.
        read = file.read()                                          # Leerlo.
        tiempo = np.double(st.unpack("d"*int(len(read)/8),read))    # Convertir a ndarray.
    with open(filevalues,"rb") as file:                             # Abrir archivo del potencial.
        read = file.read()                                          # Leerlo.
        values = np.double(st.unpack("d"*int(len(read)/8),read))    # Convertir a ndarray.
    return tiempo,values