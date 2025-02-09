# ------------------------------------------------------
# ---------------------- main.py -----------------------
# ------------------------------------------------------
from PyQt5.QtWidgets import *
from PyQt5.uic import loadUi
import numpy as np
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
from PyQt5.QtCore import *
from PyQt5.QtGui import *

import modulo as mod
import imagenes

global dict_graficadas, dict_datos, labels, colors

# Creacion de variables globales

labels = []
colors = {"Euler Forward": "b", "Euler Backward": "r", "Euler Modificado": "y", "Runge-Kutta 2": "g", "Solve_ivp": "k",
          "Runge-Kutta 4": "c", "Odeint": "m"}
dict_graficadas = {}
dict_datos = {}


# Llenar diccionarios
def crear_dicts():
    for l in ["Euler Forward", "Euler Backward", "Euler Modificado", "Runge-Kutta 2", "Runge-Kutta 4", "Odeint",
              "Solve_ivp"]:
        dict_graficadas[l] = 0
        dict_datos[l] = 0


crear_dicts()


def poderosa(fun, v, h, m, n, temp, tf, c1, cbfija, cbvari, c2, t10, t1f, t20, t2f):
    """
        La funcion poderosa se encarga de generar los datos para graficar, indenpendiente de los datos iniciales

        :param fun: Nombre del metodo numerico
        :param v: voltaje inicial
        :param h: porcentaje de puerta h
        :param m: porcentaje de puerta m
        :param n: porcentaje de puerta n
        :param temp: temperatura
        :param tf: tiempo de estimulacion total
        :param c1: corriente 1
        :param cbfija: valor del chechBox de corriente fija
        :param cbvari: valor del chechBox de corriente variable
        :param c2: corriente 2
        :param t10: tiempo inicial para corriente 1
        :param t1f: tiempo final para corriente 1
        :param t20: tiempo inicial para corriente 2
        :param t2f: tiempo final para corriente 2
        :return: tupla con las listas tiempo y valores  o -1 en caso de error.
        """

    condiciones = [v, n, m, h]  # V,N,M,H definicion de condiciones iniciales
    rev = -1  # variable que determina si hay error o se ejecuto con normalidad el codigo
    if cbfija:  # si se selecciono corriente fija
        rev = fun(condiciones, tf, temp, lambda t: float(c1))
    elif cbvari:  # si se selecciono corriente variable
        def Ivar(t):
            if t10 < t < t1f:
                return float(c1)
            elif t20 < t < t2f:
                return float(c2)
            return 0

        rev = fun(condiciones, tf, temp, Ivar)
    if rev == -1:  # En caso de error
        mensaje = QMessageBox()
        mensaje.setWindowTitle("Error")
        mensaje.setText("Elija un tipo de corriente")
        mensaje.setIcon(QMessageBox.Critical)
        mensaje.setStandardButtons(QMessageBox.Ok)
        mensaje.exec_()
        return -1
    else:
        return rev  # retorna tiempos y valores respectivos


class MatplotlibWidget(QMainWindow):  # inicializacion de ventana principal

    def __init__(self):
        QMainWindow.__init__(self)

        loadUi("bb.ui", self)  # Se carga el archivo .ui

        self.setWindowTitle("Interfaz para simulacion de transmision neuronal")  # Se pone titulo a la ventana
        self.setWindowIcon(QIcon(":imagenes/gorro.png"))  # Se pone icono

        # Asignacion de funciones a botones de metodos numericos
        self.BotonFOR.clicked.connect(self.forw)
        self.BotonBAC.clicked.connect(self.back)
        self.BotonMOD.clicked.connect(self.modi)
        self.BotonRK2.clicked.connect(self.rk2)
        self.BotonRK4.clicked.connect(self.rk4)
        self.BotonODE.clicked.connect(self.odeint)
        self.BotonIVP.clicked.connect(self.ivp)
        self.BotonLIM.clicked.connect(self.limpiar)

        # Todoo lo que tiene que ver con la configuracion del parametro H

        # inicializar valores del slicer y conectar a funcion cuando se cambia

        self.sh.setMinimum(0)
        self.sh.setMaximum(100)
        self.sh.setValue(65)
        self.sh.setTickInterval(10)
        self.sh.valueChanged.connect(self.sh_change)

        # Inicializar valores, aliniamiento y funcion al label del valor de h

        self.leh.setText(str(65))
        self.leh.setAlignment(Qt.AlignCenter)
        self.leh.editingFinished.connect(self.leh_change)
        self.lah.setAlignment(Qt.AlignCenter)

        # Todoo lo que tiene que ver con la configuracion del parametro M

        # inicializar valores del slicer y conectar a funcion cuando se cambia

        self.sm.setMinimum(0)
        self.sm.setMaximum(100)
        self.sm.setValue(5)
        self.sm.setTickInterval(10)
        self.sm.valueChanged.connect(self.sm_change)

        # Inicializar valores, aliniamiento y funcion al label del valor de m

        self.lem.setText(str(5))
        self.lem.setAlignment(Qt.AlignCenter)
        self.lem.editingFinished.connect(self.lem_change)
        self.lam.setAlignment(Qt.AlignCenter)

        # Todoo lo que tiene que ver con la configuracion del parametro N

        # inicializar valores del slicer y conectar a funcion cuando se cambia

        self.sn.setMinimum(0)
        self.sn.setMaximum(100)
        self.sn.setValue(30)
        self.sn.setTickInterval(100)
        self.sn.valueChanged.connect(self.sn_change)

        # Inicializar valores, aliniamiento y funcion al label del valor de n

        self.len.setText(str(30))
        self.len.setAlignment(Qt.AlignCenter)
        self.len.editingFinished.connect(self.len_change)
        self.lan.setAlignment(Qt.AlignCenter)

        # Aliniacion de labels de parametros
        self.c1_la.setAlignment(Qt.AlignRight)
        self.c2_la.setAlignment(Qt.AlignRight)
        self.v0_la.setAlignment(Qt.AlignRight)
        self.temp_la.setAlignment(Qt.AlignRight)

        # Agregar acciones de importar y exportar
        self.actionIm.triggered.connect(self.getBIN)
        self.actionEx.triggered.connect(self.guarBIN)

        # conectar los checkbox de corriente a funciones
        self.cBoxFIJA.stateChanged.connect(self.fija_change)
        self.cBoxVARI.stateChanged.connect(self.vari_change)

        # Carga la seccion de la grafica
        self.addToolBar(NavigationToolbar(self.grafica.canvas, self))

        # colocar imagen de fondo
        self.imagen.setPixmap(
            QPixmap(":imagenes/pulpo.png").scaled(100, 100, Qt.KeepAspectRatio, Qt.SmoothTransformation))

    def generalizar(self, fun, nombre, color):
        """
        generaliza el codigo de obtencion de datos y graficacion
        :param fun: nombre del metodo numerico
        :param nombre: nombre del metodo en el diccionario
        :param color: color para graficar
        :return: grafica en el canvas dispuesto
        """
        dict_graficadas[nombre] = 1

        v = float(self.v0_val.text())
        h = self.sh.value() / 100
        m = self.sm.value() / 100
        n = self.sn.value() / 100
        temp = float(self.temp_val.text())
        tf = float(self.tt_val.text())

        valores = poderosa(fun, v, h, m, n, temp, tf, float(self.c1.text()), self.cBoxFIJA.isChecked(),
                           self.cBoxVARI.isChecked(), float(self.c2.text()),
                           float(self.t10.text()), float(self.t1f.text()), float(self.t20.text()),
                           float(self.t2f.text()))

        if valores != -1:
            self.graficar(valores, nombre, color)
            dict_datos[nombre] = valores

    def graficar(self, valores, nombre, color):
        """
        grafica los valores, asigna titulos y legendas
        :param valores:  tupla con tiempos y valores del metodo
        :param nombre: nombre del metodo usado
        :param color: color a graficar
        :return: retorna grafica en el canvas dispuesto
        """
        if nombre not in labels:
            labels.append(nombre)
            self.grafica.canvas.axes.plot(*valores, color)
            self.grafica.canvas.axes.legend(labels, loc='upper right')
            self.grafica.canvas.axes.set_title('Estimulacion neuronal')
            self.grafica.canvas.axes.set_xlabel('Tiempo (ms)')
            self.grafica.canvas.axes.set_xlabel('Voltaje (V)')
            self.grafica.canvas.draw()

    def showEvent(self, event) -> None:

        """
        Al iniciar se mueven ciertos objetos
        """
        self.c2.move(500, 500)
        self.c2_la.move(500, 500)
        self.te_la.move(500, 500)
        self.t10.move(500, 500)
        self.t1f.move(500, 500)
        self.t20.move(500, 500)
        self.t2f.move(500, 500)
        self.ms1.move(500, 500)
        self.ms2.move(500, 500)

    def closeEvent(self, event) -> None:
        """
        Al cerrar se presenta message box para preguntar si se quiere exportar, solo aparece si hay datos graficados
        """
        for v in dict_datos.values():
            if v != 0:
                mensaje = QMessageBox()
                mensaje.setWindowTitle("Ya te vas?")
                mensaje.setText("Te gustaria exportar los datos graficados?")
                mensaje.setIcon(QMessageBox.Question)
                mensaje.setStandardButtons(QMessageBox.Yes | QMessageBox.No)

                def salida(boton):
                    if "Yes" in boton.text():
                        self.guarBIN()
                    else:
                        pass

                mensaje.buttonClicked.connect(salida)
                mensaje.exec_()

    def fija_change(self):
        """
        Se mueven objetos de corriente variable y se evita que ambos check box se activen al tiempo
        """
        if self.cBoxFIJA.isChecked() and not (self.cBoxVARI.isChecked()):
            self.c2.move(500, 500)
            self.c2_la.move(500, 500)
            self.te_la.move(500, 500)
            self.t10.move(500, 500)
            self.t1f.move(500, 500)
            self.t20.move(500, 500)
            self.t2f.move(500, 500)
            self.ms1.move(500, 500)
            self.ms2.move(500, 500)
        elif self.cBoxVARI.isChecked() and self.cBoxFIJA.isChecked():
            self.cBoxVARI.setChecked(False)
            self.cBoxFIJA.setChecked(True)
            self.c2.move(500, 500)
            self.c2_la.move(500, 500)
            self.te_la.move(500, 500)
            self.t10.move(500, 500)
            self.t1f.move(500, 500)
            self.t20.move(500, 500)
            self.t2f.move(500, 500)
            self.ms1.move(500, 500)
            self.ms2.move(500, 500)

    def vari_change(self):
        """
        Se mueven objetos de corriente variable a su posicion original y se evita que ambos check box se activen al tiempo
        """
        if self.cBoxVARI.isChecked() and not (self.cBoxFIJA.isChecked()):
            self.gridLayout_4.addWidget(self.te_la, 3, 0, 1, 1)
            self.gridLayout_4.addWidget(self.c2, 5, 7, 1, 1)
            self.gridLayout_4.addWidget(self.c2_la, 5, 6, 1, 1)
            self.gridLayout_4.addWidget(self.t10, 3, 1, 1, 1)
            self.gridLayout_4.addWidget(self.t1f, 3, 3, 1, 1)
            self.gridLayout_4.addWidget(self.t20, 5, 1, 1, 1)
            self.gridLayout_4.addWidget(self.t2f, 5, 3, 1, 1)
        elif self.cBoxFIJA.isChecked() and self.cBoxVARI.isChecked():
            self.cBoxFIJA.setChecked(False)
            self.cBoxVARI.setChecked(True)
            self.gridLayout_4.addWidget(self.te_la, 3, 0, 1, 1)
            self.gridLayout_4.addWidget(self.c2, 5, 7, 1, 1)
            self.gridLayout_4.addWidget(self.c2_la, 5, 6, 1, 1)
            self.gridLayout_4.addWidget(self.t10, 3, 1, 1, 1)
            self.gridLayout_4.addWidget(self.t1f, 3, 3, 1, 1)
            self.gridLayout_4.addWidget(self.t20, 5, 1, 1, 1)
            self.gridLayout_4.addWidget(self.t2f, 5, 3, 1, 1)

    def leh_change(self):
        """
        Si cambia el valor en la linea editable se cambia en el slicer tambien
        """
        val = int(self.leh.text())

        if val != int(self.sh.value()) and 0 <= val <= 100:
            self.sh.setValue(val)

        elif not (0 <= val <= 100):
            mensaje = QMessageBox()
            mensaje.setWindowTitle("Informacion")
            mensaje.setText("El valor debe estar entre 0 y 100")
            mensaje.setIcon(QMessageBox.Critical)
            mensaje.setStandardButtons(QMessageBox.Ok)
            self.sh_change()
            mensaje.exec_()

    def lem_change(self):
        """
        Si cambia el valor en la linea editable se cambia en el slicer tambien
        """
        val = int(self.lem.text())

        if val != int(self.sm.value()) and 0 <= val <= 100:
            self.sm.setValue(val)

        elif not (0 <= val <= 100):
            mensaje = QMessageBox()
            mensaje.setWindowTitle("Informacion")
            mensaje.setText("El valor debe estar entre 0 y 100")
            mensaje.setIcon(QMessageBox.Critical)
            mensaje.setStandardButtons(QMessageBox.Ok)
            self.sm_change()
            mensaje.exec_()

    def len_change(self):
        """
        Si cambia el valor en la linea editable se cambia en el slicer tambien
        """
        val = int(self.len.text())

        if val != int(self.sn.value()) and 0 <= val <= 100:
            self.sn.setValue(val)

        elif not (0 <= val <= 100):
            mensaje = QMessageBox()
            mensaje.setWindowTitle("Informacion")
            mensaje.setText("El valor debe estar entre 0 y 100")
            mensaje.setIcon(QMessageBox.Critical)
            mensaje.setStandardButtons(QMessageBox.Ok)
            self.sn_change()
            mensaje.exec_()

    def sh_change(self):
        """
        Si cambia el valor en el slicer se cambia en la linea editable
        """
        my_value = str(self.sh.value())
        self.leh.setText(my_value)

    def sm_change(self):
        """
        Si cambia el valor en el slicer se cambia en la linea editable
        """
        my_value = str(self.sm.value())
        self.lem.setText(my_value)

    def sn_change(self):
        """
        Si cambia el valor en el slicer se cambia en la linea editable
        """
        my_value = str(self.sn.value())
        self.len.setText(my_value)

    def forw(self):
        """
        Implementacion de Euler forward
        """
        fun = mod.EulerFor
        nombre = "Euler Forward"
        color = colors[nombre]

        self.generalizar(fun, nombre, color)

    def back(self):
        """
        Implementacion de Euler Backward
        """
        fun = mod.EulerBack
        nombre = "Euler Backward"
        color = colors[nombre]

        QApplication.setOverrideCursor(Qt.BusyCursor)
        self.generalizar(fun, nombre, color)
        QApplication.restoreOverrideCursor()

    def modi(self):
        """
        Implementacion de Euler Modificado
        """
        fun = mod.EulerMod
        nombre = "Euler Modificado"
        color = colors[nombre]

        QApplication.setOverrideCursor(Qt.BusyCursor)
        self.generalizar(fun, nombre, color)
        QApplication.restoreOverrideCursor()

    def rk2(self):
        """
        Implementacion de Runge-Kutta 2
        """
        fun = mod.RungeKutta2
        nombre = "Runge-Kutta 2"
        color = colors[nombre]

        QApplication.setOverrideCursor(Qt.BusyCursor)
        self.generalizar(fun, nombre, color)
        QApplication.restoreOverrideCursor()

    def rk4(self):
        """
        Implementacion de Runge-Kutta 4
        """
        fun = mod.RungeKutta4
        nombre = "Runge-Kutta 4"
        color = colors[nombre]

        QApplication.setOverrideCursor(Qt.BusyCursor)
        self.generalizar(fun, nombre, color)
        QApplication.restoreOverrideCursor()

    def odeint(self):
        """
        Implementacion de Odeint
        """
        fun = mod.scipy_odeint
        nombre = "Odeint"
        color = colors[nombre]

        QApplication.setOverrideCursor(Qt.BusyCursor)
        self.generalizar(fun, nombre, color)
        QApplication.restoreOverrideCursor()

    def ivp(self):
        """
        Implementacion de Solve_ivp
        """
        fun = mod.scipy_solve
        nombre = "Solve_ivp"
        color = colors[nombre]

        QApplication.setOverrideCursor(Qt.BusyCursor)
        self.generalizar(fun, nombre, color)
        QApplication.restoreOverrideCursor()

    def limpiar(self):
        self.grafica.canvas.axes.clear()
        self.grafica.canvas.draw()
        crear_dicts()
        global labels
        labels = []

    def getBIN(self):

        """
        se ejecuta importar, se muestran message box para facilidad del ususario
        """

        mensaje = QMessageBox()
        mensaje.setWindowTitle("Informacion")
        mensaje.setText("Seleccione el archivo de tiempos")
        mensaje.setIcon(QMessageBox.Information)
        mensaje.exec_()

        filetiempo, _ = QFileDialog.getOpenFileName(self, 'Open file', '/home')

        mensaje = QMessageBox()
        mensaje.setWindowTitle("Informacion")
        mensaje.setText("Seleccione el archivo de valores")
        mensaje.setIcon(QMessageBox.Information)
        mensaje.exec_()

        filevalues, _ = QFileDialog.getOpenFileName(self, 'Open file', '/home')

        if filetiempo != "" and filevalues != "":
            for i in dict_graficadas.keys():
                if i in filetiempo:
                    dict_graficadas[i] = 1
                    self.graficar(mod.importar(filetiempo, filevalues), i, colors[i])

        else:
            mensaje = QMessageBox()
            mensaje.setWindowTitle("Error")
            mensaje.setText("No selecciono uno o ninguno de los archivos")
            mensaje.setIcon(QMessageBox.Critical)
            mensaje.setStandardButtons(QMessageBox.Ok)
            mensaje.exec_()

    def guarBIN(self):

        """
        se limpia el diccionario de datos y luego se usa la funcion exportar
        """

        lista = []

        for k, v in dict_datos.items():
            if v == 0:
                lista.append(k)

        for i in lista:
            del dict_datos[i]

        mod.exportar(dict_datos)


app = QApplication([])
window = MatplotlibWidget()
window.show()
app.exec_()
