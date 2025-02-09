from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvas

import matplotlib.pyplot as plt


class grafica(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        plt.style.use('bmh')
        fig = plt.figure()

        self.canvas = FigureCanvas(fig)

        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.setLayout(vertical_layout)