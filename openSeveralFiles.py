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
        self.r_pos = []
        self.r_press = []
        self.filenames = []

        self.movement = [0]
        self.floatingVar = [0 for x in range(20)]
        self.floatVarSpd = [0 for x in range(20)]

        self.moreThanVarY = []
        self.moreThanVarX = []
        self.lat = []
        self.lon = []
        self.speed = []

        self.overMarkX = []
        self.overMarkY = []

        self.csvFileName = ''

        self.open_file()
        #self.processCalibration()
        #self.prepare_data()
        self.push_to_one_zero()
        self.show_data()
        #self.write_to_csv()


    def open_file(self):
        prevSpd = 0.0
        prevTime = 0
        self.filenames = QFileDialog.getOpenFileNames(self,'Open File', 'C:\\Users\\ADiKo\\Desktop\\', 'txt file (*.txt)')[0]
        print(self.filenames)
        #self.csvFileName = fname.split('.')[0] + ".csv"
        try:
            for fname in self.filenames:
                f = open(fname, 'r')
                temp_timings = []
                temp_l_pos = []
                temp_l_press = []
                temp_r_pos = []
                temp_r_press = []
                print("file opened")
                for line in f:
                    if line and len(line) > 5:
                        l = line.split(',')
                        if l[0] is 'a':
                            temp_timings.append(int(l[1]))
                            temp_l_pos.append(int(l[2]))
                            temp_l_press.append(int(l[3]))
                            temp_r_pos.append(int(l[4]))
                            temp_r_press.append(int(l[5]))

                f.close()
                self.timings.append(temp_timings)
                self.l_pos.append(temp_l_pos)
                self.l_press.append(temp_l_press)
                self.r_pos.append(temp_r_pos)
                self.r_press.append(temp_r_press)

        except Exception as ex:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)

        for timings_array in self.timings:
            delta = timings_array[0]
            for j in range(len(timings_array)):
                timings_array[j] -= delta


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

        q = self.l_pos.index(mass[temp][1])+1
        w = self.l_pos.index(mass[int(len(mass) - temp)][1])
        self.x2_lin = self.timings[q:w+1]
        self.y2_lin = [k1 * i + o1 for i in self.x2_lin]

        q = self.l_pos.index(mass[int(len(mass) - temp)][1])+1
        w = self.l_pos.index(mass[len(mass) - 1][1])
        self.x3_par = self.timings[q:w+1]
        self.y3_par = [a2 * i * i + b2 * i + c2 for i in self.x3_par]

        plp_approx = {'x': [], 'y': []}
        plp_approx['x'].extend(self.x1_par)
        plp_approx['x'].extend(self.x2_lin)
        plp_approx['x'].extend(self.x3_par)
        plp_approx['y'].extend(self.y1_par)
        plp_approx['y'].extend(self.y2_lin)
        plp_approx['y'].extend(self.y3_par)
        print("Len of original:", len(self.timings))

        (k1, b1) = get_linear_coeff(mass[0][0], mass[0][1], mass[temp][0], mass[temp][1])
        (k2, b2) = get_linear_coeff(mass[temp][0], mass[temp][1], mass[int(len(mass) - temp)][0], mass[int(len(mass) - temp)][1])
        (k3, b3) = get_linear_coeff(mass[int(len(mass) - temp)][0], mass[int(len(mass) - temp)][1], mass[len(mass) - 1][0], mass[len(mass) - 1][1])

        q = self.l_pos.index(mass[0][1])
        w = self.l_pos.index(mass[temp][1])
        x1 = self.timings[q:w+1]
        y1 = [k1 * i + b1 for i in x1]

        q = self.l_pos.index(mass[temp][1]) + 1
        w = self.l_pos.index(mass[int(len(mass) - temp)][1])
        x2 = self.timings[q:w+1]
        y2 = [k2 * i + b2 for i in x2]

        q = self.l_pos.index(mass[int(len(mass) - temp)][1]) + 1
        w = self.l_pos.index(mass[len(mass) - 1][1])
        x3 = self.timings[q:w+1]
        y3 = [k3 * i + b3 for i in x3]

        lll_approx = {'x': [], 'y': []}
        lll_approx['x'].extend(x1)
        lll_approx['x'].extend(x2)
        lll_approx['x'].extend(x3)
        lll_approx['y'].extend(y1)
        lll_approx['y'].extend(y2)
        lll_approx['y'].extend(y3)

        (k1, b1) = get_linear_coeff(mass[0][0], mass[0][1], mass[len(mass) - 1][0], mass[len(mass) - 1][1])
        q = self.l_pos.index(mass[0][1])
        w = self.l_pos.index(mass[len(mass) - 1][1])
        x1 = self.timings[q:w+1]
        y1 = [k1 * i + b1 for i in x1]
        l_approx = {'x': [], 'y': []}
        l_approx['x'].extend(x1)
        l_approx['y'].extend(y1)

        temp = int(len(mass) / 2)  # вообще 3, но пох)
        (a1, b1, c1) = get_parabola_coeff(mass[0][0], mass[0][1], mass[int(temp / 2)][0], mass[int(temp / 2)][1],
                                          mass[temp][0], mass[temp][1])
        (a2, b2, c2) = get_parabola_coeff(mass[temp][0], mass[temp][1], mass[int(len(mass) - temp / 2)][0], mass[int(len(mass) - temp / 2)][1],
                                          mass[len(mass) - 1][0], mass[len(mass) - 1][1])

        q = self.l_pos.index(mass[0][1])
        w = self.l_pos.index(mass[temp][1])
        x1 = self.timings[q:w + 1]
        y1 = [a1 * i * i + b1 * i + c1 for i in x1]

        q = self.l_pos.index(mass[temp][1]) + 1
        w = self.l_pos.index(mass[len(mass) - 1][1])
        x2 = self.timings[q:w + 1]
        y2 = [a2 * i * i + b2 * i + c2 for i in x2]

        pp_approx = {'x': [], 'y': []}
        pp_approx['x'].extend(x1)
        pp_approx['x'].extend(x2)
        pp_approx['y'].extend(y1)
        pp_approx['y'].extend(y2)

        temp = int(len(mass) / 3)
        (k1, b1) = get_linear_coeff(mass[0][0], mass[0][1], mass[int(len(mass) - temp)][0], mass[int(len(mass) - temp)][1])
        (a2, b2, c2) = get_parabola_coeff(mass[int(len(mass) - temp)][0], mass[int(len(mass) - temp)][1],
                                          mass[int(len(mass) - temp / 2)][0],
                                          mass[int(len(mass) - temp / 2)][1], mass[len(mass) - 1][0],
                                          mass[len(mass) - 1][1])
        q = self.l_pos.index(mass[0][1])
        w = self.l_pos.index(mass[int(len(mass) - temp)][1])
        x1 = self.timings[q:w + 1]
        y1 = [k1 * i + b1 for i in x1]

        q = self.l_pos.index(mass[int(len(mass) - temp)][1]) + 1
        w = self.l_pos.index(mass[len(mass) - 1][1])
        x2 = self.timings[q:w + 1]
        y2 = [a2 * i * i + b2 * i + c2 for i in x2]

        lp_approx = {'x': [], 'y': []}
        lp_approx['x'].extend(x1)
        lp_approx['x'].extend(x2)
        lp_approx['y'].extend(y1)
        lp_approx['y'].extend(y2)


        plt.plot(self.timings, self.l_pos, 'r')
        plt.plot(plp_approx['x'], plp_approx['y'], 'g')
        plt.plot(lll_approx['x'], lll_approx['y'], 'b')
        plt.plot(l_approx['x'], l_approx['y'], 'k')
        plt.plot(pp_approx['x'], pp_approx['y'], 'y')
        plt.plot(lp_approx['x'], lp_approx['y'], 'c')

        #print("len:", len(plp_approx['x']), len(lll_approx['x']), len(l_approx['x']), len(pp_approx['x']))

        print("err plp:", get_error_meaning(self.timings, self.l_pos, plp_approx['x'], plp_approx['y']))
        print("err lll:", get_error_meaning(self.timings, self.l_pos, lll_approx['x'], lll_approx['y']))
        print("err l:", get_error_meaning(self.timings, self.l_pos, l_approx['x'], l_approx['y']))
        print("err pp:", get_error_meaning(self.timings, self.l_pos, pp_approx['x'], pp_approx['y']))
        print("err lp:", get_error_meaning(self.timings, self.l_pos, lp_approx['x'], lp_approx['y']))
        plt.grid(True)
        plt.show()

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



    def push_to_one_zero(self):
        if len(self.l_pos) < 2:
            return
        list_of_indx = []
        list_of_indx2 = []
        for pos in self.l_pos[0]:
            got = 0
            for pos_array in self.l_pos:
                if pos in pos_array:
                    got += 1
            if got == len(self.l_pos):
                for pos_array in self.l_pos:
                    list_of_indx.append(pos_array.index(pos))
                break

        for i in range(len(self.l_pos[0]) - 1,-1,-1):
            got = 0
            for pos_array in self.l_pos:
                if self.l_pos[0][i] in pos_array:
                    got += 1
            if got == len(self.l_pos):
                for pos_array in self.l_pos:
                    list_of_indx2.append(pos_array.index(self.l_pos[0][i]))
                break

        for i in range(len(self.l_pos)):
            self.timings[i] = self.timings[i][list_of_indx[i]:list_of_indx2[i]]
            self.l_pos[i] = self.l_pos[i][list_of_indx[i]:list_of_indx2[i]]
            self.l_press[i] = self.l_press[i][list_of_indx[i]:list_of_indx2[i]]

        for i in range(len(self.l_pos)):
            delta = self.timings[i][0]
            for j in range(len(self.timings[i])):
                self.timings[i][j] -= delta

        print("indx:", list_of_indx)




    def show_data(self):
        colors = ['r', 'g', 'b', 'c', 'y', 'k']
        for i in range(len(self.timings)):
            plt.figure(1)
            plt.plot(self.timings[i], self.l_pos[i], colors[i])
            plt.figure(2)
            plt.plot(self.timings[i], self.l_press[i], colors[i])
            print(self.filenames[i], ":", colors[i], "time:", self.timings[i][-1], "mean pressure:", mean(self.l_press[i]))
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
    return a, b, c

def get_linear_coeff(x1, y1, x2, y2):
    k = 0
    b = 0
    k = (y2 - y1) / (x2 - x1)
    b = y1 - k * x1
    return k, b

def get_error_meaning(x_ref, y_ref, x, y):
    err = 0
    for i in range(len(x)):
        index = x_ref.index(x[i])
        err += (y_ref[index] - y[i]) ** 2

    err = sqrt(err)
    return err

if __name__ == '__main__':
    print(get_parabola_coeff(0, 5, 1, 6, 2, 9))
    print(get_error_meaning([1,2,3,4,5,6],[5,5,5,5,5,5],[2,3,4,5], [1,2,3,4]))
    app = QApplication(sys.argv)
    ex = DataVisualisation()
   # sys.exit(app.exec_())