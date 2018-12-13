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
        self.navinfo = []

        self.l_pos_filtered = []

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

        self.pos_cm = []

        self.csvFileName = ''

        self.open_file()
        print("Variance is", variance(self.l_pos), variance(self.r_pos), variance(self.l_press), variance(self.r_press))
        plt.figure(1)
        hist, bins = np.histogram(self.l_pos, bins=int(len(self.l_pos) ** (1/3.0)))
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        #plt.show()
        plt.figure(2)
        hist, bins = np.histogram(self.r_pos, bins=int(len(self.r_pos) ** (1 / 3.0)))
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        #plt.show()
        plt.figure(3)
        hist, bins = np.histogram(self.l_press, bins=int(len(self.l_press) ** (1 / 3.0)))
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        #plt.show()
        plt.figure(4)
        hist, bins = np.histogram(self.r_press, bins=int(len(self.r_press) ** (1 / 3.0)))
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        plt.show()

        #self.calculate_acceleration()
        #self.processCalibration()

        #self.prepare_data()
        #self.filter_position()
        #self.show_data()
        #self.write_to_csv()
    def calculate_acceleration(self):
        clearence_range = 0.13
        max_pos = 3100
        min_pos = 1900

        speed = [0]
        acceleration = [0]
        acc_2 = [0]
        pos = []
        time = []
        filtered_with_var = []
        filtered_with_var2 = []

        for p in self.l_pos:
            pos.append((p - min_pos) / (max_pos - min_pos) * clearence_range)

        for t in self.timings:
            time.append(t / 1000)

        for i in range(1, len(pos)):
            try:
                speed.append((pos[i] - pos[i-1]) / (time[i] - time[i-1]))
            except Exception as ex:
                template = "An exception of type {0} occurred. Arguments:\n{1!r}"
                message = template.format(type(ex).__name__, ex.args)
                print(message, "timings", i, time[i], time[i-1])

        for i in range(1, len(pos)):
            try:
                acceleration.append((speed[i] - speed[i-1]) / (time[i] - time[i-1]))
            except Exception as ex:
                template = "An exception of type {0} occurred. Arguments:\n{1!r}"
                message = template.format(type(ex).__name__, ex.args)
                print(message, "timings", i, time[i], time[i-1])

        for i in range(1, len(pos)):
            try:
                acc_2.append((acceleration[i] - acceleration[i-1]) / (time[i] - time[i-1]) * 10)
            except Exception as ex:
                template = "An exception of type {0} occurred. Arguments:\n{1!r}"
                message = template.format(type(ex).__name__, ex.args)
                print(message, "timings", i, time[i], time[i-1])

        fs = 40.0
        order = 5
        cutoff = 1
        filtered = butter_lowpass_filter(self.l_pos, cutoff, fs, order)

        pos_f = []

        for p in filtered:
            pos_f.append((p - min_pos) / (max_pos - min_pos) * clearence_range)

        for i in range(len(pos)-20):
            var = sqrt(variance(pos[i:i+20]))
            filtered_with_var.append(pos_f[i+20] + 2*var)
            filtered_with_var2.append(pos_f[i+20] - 2*var)

        print("len: ", len(self.timings), len(speed))
        thr = [10 for i in time]

        thr1 = [1250 for i in time]
        thr2 = [1770 for i in time]
        thr3 = [6100 for i in time]

        plt.figure(1)
        ax1 = plt.subplot(411)
        plt.plot(time, pos, 'r')
        plt.plot(time[100:], pos_f[100:], 'b')
        plt.plot(time[100:], filtered_with_var[80:], 'y')
        plt.plot(time[100:], filtered_with_var2[80:], 'y')
        #plt.setp(ax1.get_xticklabels(), fontsize=6)
        plt.ylabel("Позиция, м")
        plt.grid(True)

        ax2 = plt.subplot(412, sharex=ax1)
        plt.plot(time, speed, 'b')
        #plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel("Скорость, м/с")
        plt.grid(True)

        ax3 = plt.subplot(413, sharex=ax1)
        plt.plot(time, acceleration, 'g')
        plt.plot(time, thr, 'y')
        #plt.setp(ax3.get_xticklabels(), visible=False)
        plt.ylabel("Ускорение, м/с^2")
        plt.grid(True)

        ax4 = plt.subplot(414, sharex=ax1)
        plt.plot(time, acc_2, 'k')
        #plt.plot(time, thr1, 'y')
        #plt.plot(time, thr2, 'y')
        plt.plot(time, thr3, 'y')
        #plt.setp(ax3.get_xticklabels(), visible=False)
        plt.ylabel("Рывок, cм/с^3")
        plt.xlabel("Время, с")
        plt.grid(True)
        plt.show()

    def filter_position(self):
        dang_x = []
        dang_y = []

        fs = 40.0
        order = 5
        cutoff = 1
        filtered = []
        filtered = butter_lowpass_filter(self.l_pos, cutoff, fs, order)
        #filtered_mov = butter_lowpass_filter(self.movement, 2.0, fs, order)
        print("len 1:", len(self.l_pos), "len 2:", len(filtered))

        triggered = False
        start_ind = 0
        end_ind = 0
        start_pos = 0
        p_time = 0
        start_filtered = 0
        road_imp = {'start_ind': 0,
                    'start_pos': 0,
                    'max_ind': 0,
                    'max_pos': 0,
                    'max_pos_time': 0,
                    'end_ind': 0,
                    'navinfo': (0,0)
                    }
        road_imp_list = []
        for i in range(len(self.l_pos)):
            temp = abs(self.l_pos[i] - filtered[i])
            if temp > 200:
                if not triggered:
                    triggered = True
                    road_imp['start_ind'] = i
                    road_imp['start_pos'] = self.l_pos[i]
                    start_filtered = filtered[i]
            else:
                if triggered:
                    triggered = False
                    road_imp['end_ind'] = i
                    road_imp['max_pos'] = road_imp['start_pos']
                    road_imp['max_ind'] = road_imp['start_ind']
                    road_imp['max_pos_time'] = self.timings[road_imp['start_ind']]
                    road_imp['navinfo'] = self.navinfo[road_imp['start_ind']]

                    for h in range(road_imp['start_ind'], road_imp['end_ind']+1, 1):
                        if abs(filtered[h] - self.l_pos[h]) > abs(filtered[h] - road_imp['max_pos']):
                            road_imp['max_pos'] = self.l_pos[h]
                            road_imp['max_ind'] = h
                            road_imp['max_pos_time'] = self.timings[h]
                            road_imp['navinfo'] = self.navinfo[h]
                    road_imp_list.append(road_imp.copy())

        for p in road_imp_list:
            print(p,)
        plt.plot(self.timings, self.l_pos, 'r')
        plt.plot(self.timings[100:], filtered[100:], 'b')
        dang_x = [o['max_pos_time'] for o in road_imp_list]
        dang_y = [o['max_pos'] for o in road_imp_list]
        plt.scatter(dang_x, dang_y, color='k')
        #plt.plot(self.timings, self.movement, 'g')
        #plt.plot(self.timings[100:], filtered_mov[100:], 'k')
        self.write_to_csv(road_imp_list)
        plt.grid(True)
        plt.show()

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
                                self.navinfo.append(last_navinfo)
                        else:
                            self.timings.append(int(l[1]))
                            self.l_pos.append(int(l[2]))
                            self.l_press.append(int(l[3]))
                            self.r_pos.append(int(l[4]))
                            self.r_press.append(int(l[5]))
                            self.navinfo.append(last_navinfo)

                    if l[0] is 'b':
                        lat1 = float(l[1][:2])
                        lat2 = float(l[1][2:]) * 100 / 60 / 100;
                        lat = lat1 + lat2;
                        if l[2] is 'S':
                            lat = -lat

                        lon1 = float(l[3][:2])
                        lon2 = float(l[3][2:]) * 100 / 60 / 100;
                        lon = lon1 + lon2;
                        if l[4] is 'W':
                            lon = -lon

                        self.lat.append(lat)
                        self.lon.append(lon)
                        self.speed.append(float(l[5]))
                        last_navinfo = (lat, lon)


            print("Number of GPS markers:", len(self.lat))
            for i in range(5):
                print(self.lat[i], self.lon[i], self.speed[i])
            f.close()

        except Exception as ex:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)

        delta = self.timings[0]
        print("D:", delta)
        for i in range(len(self.timings)):
            self.timings[i] -= delta
        print("done", self.timings[26595:])

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

        print("coeffs:", (a1, b1, c1), (a2, b2, c2))
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


        plt.figure(1)
        plt.plot(self.timings, self.l_pos, 'k')
        plt.plot(plp_approx['x'], plp_approx['y'], 'k--')
        plt.title('График изменения клиренса')
        plt.legend(('Эксперимент', 'Аппроксимация'))
        plt.xlabel("Время, мс")
        plt.ylabel("Позиция")
        plt.grid(True)

        plt.figure(2)
        plt.plot(self.timings, self.l_pos, 'k')
        plt.plot(lll_approx['x'], lll_approx['y'], 'k--')
        plt.title('График изменения клиренса')
        plt.legend(('Эксперимент', 'Аппроксимация'))
        plt.xlabel("Время, мс")
        plt.ylabel("Позиция")
        plt.grid(True)

        plt.figure(3)
        plt.plot(self.timings, self.l_pos, 'k')
        plt.plot(l_approx['x'], l_approx['y'], 'k--')
        plt.title('График изменения клиренса')
        plt.legend(('Эксперимент', 'Аппроксимация'))
        plt.xlabel("Время, мс")
        plt.ylabel("Позиция")
        plt.grid(True)

        plt.figure(4)
        plt.plot(self.timings, self.l_pos, 'k')
        plt.plot(pp_approx['x'], pp_approx['y'], 'k--')
        plt.title('График изменения клиренса')
        plt.legend(('Эксперимент', 'Аппроксимация'))
        plt.xlabel("Время, мс")
        plt.ylabel("Позиция")
        plt.grid(True)

        plt.figure(5)
        plt.plot(self.timings, self.l_pos, 'k')
        plt.plot(lp_approx['x'], lp_approx['y'], 'k--')
        plt.title('График изменения клиренса')
        plt.legend(('Эксперимент', 'Аппроксимация'))
        plt.xlabel("Время, мс")
        plt.ylabel("Позиция")
        plt.grid(True)

        #print("len:", len(plp_approx['x']), len(lll_approx['x']), len(l_approx['x']), len(pp_approx['x']))

        print("err plp:", get_error_meaning(self.timings, self.l_pos, plp_approx['x'], plp_approx['y']))
        print("err lll:", get_error_meaning(self.timings, self.l_pos, lll_approx['x'], lll_approx['y']))
        print("err l:", get_error_meaning(self.timings, self.l_pos, l_approx['x'], l_approx['y']))
        print("err pp:", get_error_meaning(self.timings, self.l_pos, pp_approx['x'], pp_approx['y']))
        print("err lp:", get_error_meaning(self.timings, self.l_pos, lp_approx['x'], lp_approx['y']))
        plt.grid(True)
        plt.xlabel("Время, мс")
        plt.ylabel("Позиция")
        plt.show()



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
        plt.plot(self.timings, self.l_pos, 'k')
        plt.title("График уменьшения клиренса")
        plt.xlabel("Время, мс")
        plt.ylabel("Позиция")
        #plt.plot(self.timings, self.l_press, 'b')
        #plt.plot(self.timings, self.movement, 'g')
        #plt.plot(self.timings, self.floatingVar, 'b')
        #plt.plot(self.timings, self.floatVarSpd, 'k')
        #plt.scatter(self.overMarkX, self.overMarkY, c='k')
        #plt.scatter(self.moreThanVarX, self.moreThanVarY, c='c')

        #plt.plot(x, y1, 'y')
        #plt.plot(x, y2, 'y')

        #plt.plot(self.x1_par, self.y1_par, 'g')
        #plt.plot(self.x2_lin, self.y2_lin, 'k')
        #plt.plot(self.x3_par, self.y3_par, 'b')

        plt.grid(True)
        plt.show()

    def write_to_csv(self, data):
        try:
            #print(self.csvFileName)
            f = open(self.csvFileName, 'w')
            f.write('Coordinates\n')

            for imp in data:
                if imp['navinfo'][0] != 0.0:
                    f.write('\"' + str(imp['navinfo'][0]) + ',' + str(imp['navinfo'][1]) + '\"\n')
            f.close()
        except Exception as ex:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)



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

    err = sqrt(err) / len(x)
    return err

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    print("b coeff:", b)
    print("a coeff:", a)
    y = lfilter(b, a, data)
    return y

if __name__ == '__main__':
    print(get_parabola_coeff(0, 5, 1, 6, 2, 9))
    print(get_error_meaning([1,2,3,4,5,6],[5,5,5,5,5,5],[2,3,4,5], [1,2,3,4]))
    app = QApplication(sys.argv)
    ex = DataVisualisation()
    sys.exit(app.exec_())