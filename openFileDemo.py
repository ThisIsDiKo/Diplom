from PyQt5.QtWidgets import QFileDialog, QMainWindow, QApplication
import sys
import matplotlib.pyplot as plt
from statistics import median, variance, mean
from math import sqrt


class DataVisualisation(QMainWindow):
    def __init__(self):
        super().__init__()

        self.l_pos = []
        self.timings = []
        self.l_press = []
        self.movement = [0]
        self.floatingVar = [0 for x in range(20)]
        self.floatVarSpd = [0 for x in range(20)]

        self.moreThanVarY = []
        self.moreThanVarX = []
        self.lat = []
        self.lon = []
        self.markersSpd = []

        self.overMarkX = []
        self.overMarkY = []

        self.csvFileName = ''

        self.open_file()
        self.processCalibration()
        #self.prepare_data()
        self.show_data()
        #self.write_to_csv()


    def open_file(self):
        prevSpd = 0.0
        prevTime = 0
        fname = QFileDialog.getOpenFileName(self,'Open File', 'C:\\Users\\ADiKo\\Desktop\\', 'txt file (*.txt)')[0]
        print(fname)
        self.csvFileName = fname.split('.')[0] + ".csv"
        try:
            f = open(fname, 'r')
            print("file opened")
            for line in f:
                print(line, end='')
                l = line.split(',')
                if l[0] is 'a':
                    self.timings.append(int(l[1]))
                    self.l_pos.append(int(l[2]))
                    self.l_press.append(int(l[3]))

                    prevSpd = float(l[6])
                    #prevTime = int(l[1])
                if l[0] is 'b':
                    lat1 = float(l[1].split(':')[1][:2])
                    lat2 = float(l[1].split(':')[1][2:]) * 100 / 60 / 100;
                    lat = lat1 + lat2;
                    if l[2] is 'S':
                        lat = -lat

                    lon1 = float(l[3].split(':')[1][:2])
                    lon2 = float(l[3].split(':')[1][2:]) * 100 / 60 / 100;
                    lon = lon1 + lon2;
                    if l[4] is 'W':
                        lon = -lon
                    self.lat.append(lat)
                    self.lon.append(lon)
                    self.markersSpd.append(prevSpd)
            print("NUmber of GPS markers:", len(self.lat))
            for i in range(5):
                print(self.lat[i], self.lon[i], self.markersSpd[i])
            f.close()
        except:
            print("Wrong file")
       # self.read_file_line_by_line(fname)

    def processCalibration(self):
        numOfIntervals = 64
        #threshold = int(4096 / numOfIntervals)
        threshold = 50
        mass = []
        #mass.append((self.timings[0], self.l_pos[0]))
        prevP = self.l_pos[0]
        for p in self.l_pos:
            if abs(prevP - p) > threshold:
                mass.append((self.timings[self.l_pos.index(p)], p))
                prevP = p

        print("len of mass:", len(mass))
        temp = int(len(mass) / 3) #вообще 3, но пох)

        (a1, b1, c1) = get_parabola_coeff(mass[0][0], mass[0][1], mass[int(temp / 2)][0], mass[int(temp / 2)][1],
                                          mass[temp][0], mass[temp][1])

        (k1, o1) = get_linear_coeff(mass[temp][0], mass[temp][1],
                                    mass[int(len(mass) - temp)][0], mass[int(len(mass) - temp)][1])

        (a2, b2, c2) = get_parabola_coeff(mass[int(len(mass) - temp)][0], mass[int(len(mass) - temp)][1],
                                          mass[int(len(mass) - temp / 2)][0],
                                          mass[int(len(mass) - temp / 2)][1], mass[len(mass) - 1][0],
                                          mass[len(mass) - 1][1])

        q = self.l_pos.index(mass[0][1])
        w = self.l_pos.index(mass[temp][1])
        self.x1_par = self.timings[q:w+1]
        self.y1_par = [a1 * i * i + b1 * i + c1 for i in self.x1_par]

        q = self.l_pos.index(mass[temp][1])
        w = self.l_pos.index(mass[int(len(mass) - temp)][1])
        self.x2_lin = self.timings[q:w+1]
        self.y2_lin = [k1 * i + o1 for i in self.x2_lin]

        q = self.l_pos.index(mass[int(len(mass) - temp)][1])
        w = self.l_pos.index(mass[len(mass) - 1][1])
        self.x3_par = self.timings[q:w+1]
        self.y3_par = [a2 * i * i + b2 * i + c2 for i in self.x3_par]

        for k in mass:
            print(k)


    def prepare_data(self):
        d = self.timings[0]
        for i in range(len(self.timings)):
            self.timings[i] -= d

        for i in range(1, len(self.l_pos)):
            self.movement.append(self.l_pos[i] - self.l_pos[i-1])

        for i in range(len(self.l_pos) - 20):
            self.floatingVar.append(sqrt(int(variance(self.l_pos[i:i+20]))))

        startSKO = sqrt(int(variance(self.movement[0:20])))
        mass = []
        for i in range(20, len(self.movement)):
            if self.movement[i] > 10 * startSKO:
                self.overMarkX.append(self.timings[i])
                self.overMarkY.append(self.movement[i])
            else:
                mass.append(self.movement[i])
                if len(mass) == 20:
                    startSKO = sqrt(int(variance(mass)))
                    mass = []
        #print("markers:", self.overMark[:20])
        for i in range(20, len(self.movement)):
            self.floatVarSpd.append(sqrt(int(variance(self.movement[i-20:i]))))
        #med = median(self.l_pos)
        med = 2727
        for i in range(len(self.l_pos)):
            self.l_pos[i] -= med

        startSKO = sqrt(int(variance(self.l_pos[0:20])))
        meanMass = int(mean(self.l_pos[0:20]))
        mass = []
        startedImp = False;
        for i in range(20, len(self.l_pos)):
            if abs(meanMass - self.l_pos[i]) > 2 * startSKO:
                if abs(self.movement[i]) > 5 * self.floatVarSpd[i]:
                    startedImp = True
                    self.moreThanVarY.append(self.l_pos[i])
                    self.moreThanVarX.append(self.timings[i])
            else:
                if startedImp:
                    startedImp = False
                    self.moreThanVarY.append(self.l_pos[i])
                    self.moreThanVarX.append(self.timings[i])

            startSKO = sqrt(int(variance(self.l_pos[i-19:i+1])))
            meanMass = mean(self.l_pos[i-19:i+1])




    def show_data(self):
        x = self.timings[:]
        y1 = [130 for i in x]
        y2 = [-130 for i in x]
        print("timings: ", self.timings[:10])
        print(len(self.timings), len(self.floatingVar))
        plt.plot(self.timings, self.l_pos, 'r')
        #plt.plot(self.timings, self.l_press, 'b')
        #plt.plot(self.timings, self.movement, 'g')
        #plt.plot(self.timings, self.floatingVar, 'b')
        #plt.plot(self.timings, self.floatVarSpd, 'k')
        #plt.scatter(self.overMarkX, self.overMarkY, c='k')
        #plt.scatter(self.moreThanVarX, self.moreThanVarY, c='c')

        #plt.plot(x, y1, 'y')
        #plt.plot(x, y2, 'y')

        plt.plot(self.x1_par, self.y1_par, 'g')
        plt.plot(self.x2_lin, self.y2_lin, 'k')
        plt.plot(self.x3_par, self.y3_par, 'b')
        plt.grid(True)
        plt.show()

    def write_to_csv(self):
        try:
            #print(self.csvFileName)
            print(len(self.lat),len(self.lon), len(self.markersSpd))
            f = open(self.csvFileName, 'w')
            f.write('Coordinates,speed\n')

            for i in range(len(self.lat)):
                f.write('\"' + str(self.lat[i]) + ',' + str(self.lon[i]) + '\",\"' + str(self.markersSpd[i])+'\"\n')
            f.close()
        except:
            print("error")




def get_parabola_coeff(x1, y1, x2, y2, x3, y3):
    a = 0
    b = 0
    c = 0

    a = (y3 - y1 - (y2 - y1) / (x2 - x1) * (x3 - x1)) / (x3**2 - x1**2 - (x1 + x2) * (x3 - x1))
    b = (y2 - y1 + a * x1 * x1 - a * x2 * x2) / (x2 - x1)
    c = y1 - a * x1 * x1 - b * x1

    return (a, b, c)

def get_linear_coeff(x1, y1, x2, y2):
    k = 0
    b = 0

    k = (y2 - y1) / (x2 - x1)
    b = y1 - k * x1

    return (k, b)

if __name__ == '__main__':
    print(get_parabola_coeff(0, 5, 1, 6, 2, 9))
    app = QApplication(sys.argv)
    ex = DataVisualisation()
   # sys.exit(app.exec_())