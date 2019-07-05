from PyQt5.QtWidgets import QFileDialog, QMainWindow, QApplication
import sys
import matplotlib.pyplot as plt
from statistics import median, variance, mean
from math import sqrt
import numpy as np
from scipy.signal import butter, lfilter, freqz


class DataVisualisation(QMainWindow):
    def __init__(self):
        super().__init__()

        self.l_pos = []
        self.timings = []
        self.l_press = []
        self.r_pos = []
        self.r_press = []


        self.working_list = []
        self.working_time = []

        self.var = []
        self.t_var = []

        self.open_file()
        self.simulate_contrast_algorithm()

        self.show_data()



    def open_file(self):
        prevSpd = 0.0
        prevTime = 0
        last_navinfo = (0,0)
        fname = QFileDialog.getOpenFileName(self,'Open File', 'C:\\Users\\ADiKo\\Desktop\\', 'txt file (*.txt)')[0]
        print(fname)
        self.csvFileName = fname.split('.')[0] + ".csv"
        try:
            f = open(fname, 'r')
            print("file opened")
            for line in f:
                print(line, end='')
                if line:
                    l = line.split(',')
                    if l[0] is 'a':
                        if self.timings:
                            if int(l[1]) - self.timings[-1] >= 25:
                                self.timings.append(int(l[1]))
                                self.l_pos.append(int(l[2]))
                                self.l_press.append(int(l[3]))
                                self.r_pos.append(int(l[4]))
                                self.r_press.append(int(l[5]))
                        else:
                            self.timings.append(int(l[1]))
                            self.l_pos.append(int(l[2]))
                            self.l_press.append(int(l[3]))
                            self.r_pos.append(int(l[4]))
                            self.r_press.append(int(l[5]))

            f.close()

        except Exception as ex:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)

        delta = self.timings[0]

        for i in range(len(self.timings)):
            self.timings[i] -= delta

        self.working_list = self.l_pos[self.timings.index(626015):self.timings.index(799455)]
        self.working_time = self.timings[self.timings.index(626015):self.timings.index(799455)]

        self.filtered = [self.working_list[0] for i in range(10)]
        self.filtered.extend(fir_lowpass_filter(self.working_time, self.working_list)[1])

        delta = self.working_time[0]

        for i in range(len(self.working_time)):
            self.working_time[i] -= delta

    def simulate_contrast_algorithm(self):
        self.var = []
        self.t_var = []

        for i in range(40, len(self.working_list), 40):
            #print(self.working_list[i-40:i])
            self.t_var.append(self.working_time[i])
            self.var.append(self.calculate_sd(self.working_list[i-40:i]))


        for i in range(len(self.var)):
            if self.var[i] < 2:
                self.var[i] = 1
            else:
                self.var[i] = 0


    def calculate_sd(self, l):
        meanL = 0
        varL = 0
        for i in range(len(l)):
            meanL += l[i]
        meanL = int(meanL / len(l))

        for i in range(len(l)):
            varL += (l[i] - meanL)**2

        varL = int(varL / (len(l) - 1))
        varL = int(sqrt(varL))
        return varL


    def show_data(self):
        plt.figure(1)
        plt.plot(self.working_time, self.working_list, 'k')
        plt.plot(self.working_time, self.filtered, 'k--')
        plt.title("График изменения высоты шасси над землей при движении автомобиля", fontsize=18)
        plt.xlabel("Время, мс")
        plt.ylabel("Позиция")
        plt.xlim(0, max(self.working_time))
        plt.legend(('Исходный сигнал', 'Сигнал на выходе ФНЧ'))
        plt.grid(True)

        plt.figure(2)
        ax1 = plt.subplot(211)
        plt.title("График изменения высоты шасси над землей при движении автомобиля", fontsize=18)
        plt.plot(self.working_time, self.working_list, 'k')
        plt.grid(True)
        plt.ylabel("Позиция")

        ax2 = plt.subplot(212, sharex=ax1)
        plt.title("Принятие решения о возможности измерения клиренса")
        plt.plot(self.t_var, self.var, 'k')
        plt.ylabel("Решение")

        plt.xlabel("Время, мс")

        plt.xlim(0, max(self.working_time))

        plt.grid(True)




        plt.show()


def fir_lowpass_filter(time, data):
    filtered_data = []
    time_filtered_data = []

    FIRCoef = [
        24691,
        26057,
        27088,
        27733,
        27952,
        27733,
        27088,
        26057,
        24691,
        23048
    ]
    DCgain = 262144

    for i in range(10, len(data)):
        y = 0
        for j in range(len(FIRCoef)):
            y += FIRCoef[j] * data[i-10+j]
        y = int(y / DCgain)

        filtered_data.append(y)
        time_filtered_data.append(time[i])
    return  time_filtered_data, filtered_data

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = DataVisualisation()
    #sys.exit(app.exec_())