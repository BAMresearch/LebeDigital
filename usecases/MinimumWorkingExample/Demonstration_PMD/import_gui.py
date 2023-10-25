import sys
import os


from PyQt5.QtGui import QFont, QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QLabel, QVBoxLayout, QPushButton, QComboBox, QWidget
from PyQt5.QtCore import pyqtSlot



class App(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Converter')
        self.setMinimumSize(400, 300)
        self.dataDict = {}

        mainWidget = QWidget(self)
        self.setCentralWidget(mainWidget)

        layout = QVBoxLayout()
        mainWidget.setLayout(layout)

        self.label1 = QLabel('Tool zum Einfügen von Daten von Excel Tabellen', self)
        font = QFont()
        font.setPointSize(14)
        font.setBold(True)
        self.label1.setFont(font)
        layout.addWidget(self.label1)
        layout.addStretch()

        self.label2 = QLabel('Welchen Datensatz wollen Sie hochladen?', self)
        layout.addWidget(self.label2)

        self.comboBox = QComboBox(self)
        self.comboBox.addItem("Mixture Design")
        self.comboBox.addItem("E-Module")
        self.comboBox.addItem("Compressive")
        layout.addWidget(self.comboBox)

        self.btn = QPushButton('Datei auswählen', self)
        self.btn.clicked.connect(self.openFileNameDialog)
        layout.addWidget(self.btn)

        self.filePathLabel = QLabel("", self)
        layout.addWidget(self.filePathLabel)

        self.doneButton = QPushButton('Fertig', self)
        self.doneButton.clicked.connect(self.closeWindow)
        layout.addStretch()
        layout.addWidget(self.doneButton)

    @pyqtSlot()
    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "Excel Files (*.xlsx *.csv *.dat)", options=options)
        if fileName:
            dataset = self.comboBox.currentText()
            self.dataDict[dataset] = fileName
            self.filePathLabel.setText(self.filePathLabel.text() + "\n" + os.path.basename(fileName))
            #self.comboBox.removeItem(self.comboBox.currentIndex())

    @pyqtSlot()
    def closeWindow(self):
        self.close()


def main():
    app = QApplication(sys.argv)
    ex = App()
    ex.show()
    app.exec_()

    # Specify the key for the conditional check
    key_to_check = 'E-Module'

    # Check if the key matches the desired condition
    if key_to_check in ex.dataDict and key_to_check == 'E-Module':
        # Expression to execute if the condition is true
        ex.dataDict[key_to_check] = ex.dataDict[key_to_check].rpartition('/')[0]



    return ex.dataDict


if __name__ == '__main__':
    data = main()
    print(data)