import sys


from PyQt5.QtGui import QFont, QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QLabel, QVBoxLayout, QPushButton, QComboBox, QWidget
from PyQt5.QtCore import pyqtSlot
import os


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
        #self.comboBox.addItem("E-Module")
        #self.comboBox.addItem("Compressive Strength")
        #self.comboBox.addItem("Mixture Design")
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

        # Create the item model
        self.item_model = QStandardItemModel(self)
        self.comboBox.setModel(self.item_model)

        # Add main items to the QComboBox
        self.module_item = QStandardItem("E-Module")
        self.mixture_item = QStandardItem("Mixture Design")
        self.compress_item = QStandardItem("Compressive Strength")

        self.item_model.appendRow(self.module_item)
        self.item_model.appendRow(self.mixture_item)
        self.item_model.appendRow(self.compress_item)

        # Add subitems only for "E-Module" and make them checkable
        mix_subitem = QStandardItem("mix")
        mix_subitem.setCheckable(True)
        specimen_subitem = QStandardItem("Specimen")
        specimen_subitem.setCheckable(True)
        self.module_item.appendRow(mix_subitem)
        self.module_item.appendRow(specimen_subitem)

    @pyqtSlot()
    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        fileName, _ = QFileDialog.getOpenFileName(self, "Datei auswählen", "",
                                                  "Excel Files (*.xlsx *.csv *.dat)", options=options)
        if fileName:
            dataset = self.comboBox.currentText()
            #selected_subitems = [self.module_item.child(row).text() for row in range(self.module_item.rowCount())
                                 #if self.module_item.child(row).checkState()]

            if dataset not in self.dataDict:
                self.dataDict[dataset] = []
                #print('check value')
                #print(selected_subitems)
            #self.dataDict[dataset].append((fileName, fileName.split("/")[-1]))
            self.dataDict[dataset].append(fileName)

            #self.dataDict[dataset] = fileName
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


    # Filter out entries with empty lists from dataDict
    #filtered_data_dict = {key: value for key, value in ex.dataDict.items() if value}

    #return filtered_data_dict

    return ex.dataDict
    #new_dict = {}
    #for key, value in ex.dataDict.items():
        #new_dict.update({key: (str(value[0]), str(value[1]))})

    #return new_dict


if __name__ == '__main__':
    data = main()
    print(data)