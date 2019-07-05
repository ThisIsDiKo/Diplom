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

        self.pos_info = {"l_pos": 0,
                         "r_pos": 0,
                         "time": 0,
                         "lat": 0,
                         "lon": 0,
                         "car_spd": 0}

        self.process_raw_info = {"l_pos":[],
                                "r_pos":[],
                                "time":[],
                                "lat":[],
                                "lon":[],
                                "car_spd": []}

        self.process_info = {"l_pos": [],
                                 "r_pos": [],
                                 "time": [],
                                 "lat": [],
                                 "lon": [],
                                 "car_spd": []}


        self.open_file()
        self.show_graphs()

        #self.prepare_data()
        #self.filter_position()
        #self.show_data()
        #self.write_to_csv()

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
        last_lat = 0
        last_lon = 0
        last_spd = 0

        fname = QFileDialog.getOpenFileName(self,'Open File', 'C:\\Users\\ADiKo\\Desktop\\', 'txt file (*.txt)')[0]
        print(fname)
        self.csvFileName = fname.split('.')[0] + ".csv"
        try:
            f = open(fname, 'r')
            print("file opened")
            for line in f:
                if line:
                    l = line.split(',')
                    #print(l)
                    if l[0] is 'a':
                        if self.process_raw_info['time']:
                            if int(l[1]) - self.process_raw_info['time'][-1] >= 25:
                                self.pos_info['time'] = int(l[1])
                                self.pos_info['l_pos'] = int(l[2])
                                self.pos_info['r_pos'] = int(l[4])

                                self.pos_info['lat'] = last_lat
                                self.pos_info['lon'] = last_lon
                                self.pos_info['car_spd'] = last_spd

                                self.process_raw_info['time'].append(self.pos_info['time'])
                                self.process_raw_info['l_pos'].append(self.pos_info['l_pos'])
                                self.process_raw_info['r_pos'].append(self.pos_info['r_pos'])

                                self.process_raw_info['lat'].append(self.pos_info['lat'])
                                self.process_raw_info['lon'].append(self.pos_info['lon'])
                                self.process_raw_info['car_spd'].append(self.pos_info['car_spd'])
                        else:
                            self.pos_info['time'] = int(l[1])
                            self.pos_info['l_pos'] = int(l[2])
                            self.pos_info['r_pos'] = int(l[4])

                            self.pos_info['lat'] = last_lat
                            self.pos_info['lon'] = last_lon
                            self.pos_info['car_spd'] = last_spd

                            self.process_raw_info['time'].append(self.pos_info['time'])
                            self.process_raw_info['l_pos'].append(self.pos_info['l_pos'])
                            self.process_raw_info['r_pos'].append(self.pos_info['r_pos'])

                            self.process_raw_info['lat'].append(self.pos_info['lat'])
                            self.process_raw_info['lon'].append(self.pos_info['lon'])
                            self.process_raw_info['car_spd'].append(self.pos_info['car_spd'])

                    #print(self.pos_info)


                    if l[0] is 'b':
                        lat1 = float(l[1][:2])
                        lat2 = float(l[1][2:]) * 100 / 60 / 100
                        lat = lat1 + lat2
                        if l[2] is 'S':
                            lat = -lat

                        lon1 = float(l[3][:2])
                        lon2 = float(l[3][2:]) * 100 / 60 / 100
                        lon = lon1 + lon2
                        if l[4] is 'W':
                            lon = -lon

                        spd = float(l[5])

                        last_lat = lat
                        last_lon = lon
                        last_spd = spd
            f.close()

        except Exception as ex:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)

        delta = self.process_raw_info['time'][0]
        print("D:", delta)
        for i in range(len(self.process_raw_info['time'])):
            self.process_raw_info['time'][i] -= delta

        print(self.process_raw_info['time'][:30])

    def adc_to_m(self, pos, clearence_range=0.13, max_pos=3100, min_pos=1900):
        pos_m = []

        for p in pos:
            pos_m.append((p - min_pos) / (max_pos - min_pos) * clearence_range)

        return pos_m

    def calculate_filtered_pos(self, pos):
        pass

    def time_ms_to_s(self, time):
        t = []
        for i in time:
            t.append(i / 1000)
        return t

    def calculate_spd(self, time, pos):
        spd = []
        spd.append(0.0)
        for i in range(len(pos) - 1):
            spd.append((pos[i+1] - pos[i]) / (time[i+1] - time[i]))
        return spd

    def calculate_accel(self, time, spd):
        accel = []
        accel.append(0.0)
        for i in range(len(spd) - 1):
            accel.append((spd[i+1] - spd[i]) / (time[i+1] - time[i]))
        return accel

    def calculate_spurt(self, time, accel):
        spurt = []
        spurt.append(0.0)
        for i in range(len(accel) - 1):
            spurt.append((accel[i+1] - accel[i]) / (time[i+1] - time[i]))
        return spurt

    def spd_to_kmh(self, spd):
        speed = []
        for s in spd:
            speed.append(1.61 * s)
        return speed

    def search_by_accel(self, time, accel, lat, lon, car_spd):
        last_imp_time = 0
        tL = []
        dang_imp = []
        print("search:", len(time), len(car_spd))
        for i in range(len(accel)):
            if abs(accel[i]) > 10:
                if time[i] - last_imp_time > 0.5:
                    tL = []
                    print(i, time[i],lat[i],lon[i], accel[i], car_spd[i])

                    tL.append(accel[i])
                    tL.append(lat[i])
                    tL.append(lon[i])
                    tL.append(car_spd[i])
                    dang_imp.append(tL)

                    last_imp_time = time[i]

        return dang_imp

    def search_by_spd(self, time, spd, lat, lon, car_spd):
        last_imp_time = 0
        tL = []
        dang_imp = []

        for i in range(len(spd)):
            if abs(spd[i]) > 0.5:
                if time[i] - last_imp_time > 0.5:
                    tL = []

                    tL.append(spd[i])
                    tL.append(lat[i])
                    tL.append(lon[i])
                    tL.append(car_spd[i])
                    dang_imp.append(tL)

                    last_imp_time = time[i]

        return dang_imp

    def search_by_spurt(self, time, spurt, lat, lon, car_spd):
        last_imp_time = 0
        tL = []
        dang_imp = []

        for i in range(len(spurt)):
            if abs(spurt[i]) > 61:
                if time[i] - last_imp_time > 0.5:
                    tL = []

                    tL.append(spurt[i])
                    tL.append(lat[i])
                    tL.append(lon[i])
                    tL.append(car_spd[i])
                    dang_imp.append(tL)

                    last_imp_time = time[i]

        return dang_imp

    def search_by_trend(self, time, data, lat, lon, car_spd):
        filtered_data = self.filter_data(data, cutoff=2)
        last_imp_time = 0
        dang_imp = []

        for i in range(40, len(data)):
            var = variance(data[i-40:i])
            var = sqrt(var)
            threshold = 4 * var
            #print(filtered_data[i],  data[i], threshold)
            if abs(filtered_data[i] - data[i]) > threshold:
                if time[i] - last_imp_time > 0.5:
                    tL = []

                    tL.append(data[i])
                    tL.append(lat[i])
                    tL.append(lon[i])
                    tL.append(car_spd[i])
                    dang_imp.append(tL)

                    last_imp_time = time[i]

        return dang_imp


    def show_graphs(self):
        self.process_info['time'] = self.time_ms_to_s(self.process_raw_info['time'])
        self.process_info['l_pos'] = self.adc_to_m(self.process_raw_info['l_pos'])
        self.process_info['l_spd'] = self.calculate_spd(self.process_info['time'], self.process_info['l_pos'])
        self.process_info['l_accel'] = self.calculate_accel(self.process_info['time'], self.process_info['l_spd'])
        #self.process_info['car_spd'] = self.spd_to_kmh(self.process_raw_info['car_spd'])
        # Т.к. данные уже в км/ч
        self.process_info['car_spd'] = self.process_raw_info['car_spd']
        self.process_info['lat'] = self.process_raw_info['lat'][:]
        self.process_info['lon'] = self.process_raw_info['lon'][:]
        self.process_info['l_spurt'] = self.calculate_spurt(self.process_info['time'], self.process_info['l_accel'])

        print(len(self.process_info['time']),
              len(self.process_info['car_spd']),
              len(self.process_info['lat']),
              len(self.process_info['lon']))

        imps_accel = self.search_by_accel(self.process_info['time'],
                                 self.process_info['l_accel'],
                                 self.process_info['lat'],
                                 self.process_info['lon'],
                                 self.process_info['car_spd'])

        imps_spd = self.search_by_spd(self.process_info['time'],
                                 self.process_info['l_spd'],
                                 self.process_info['lat'],
                                 self.process_info['lon'],
                                 self.process_info['car_spd'])

        imps_spurt = self.search_by_spurt(self.process_info['time'],
                                      self.process_info['l_spurt'],
                                      self.process_info['lat'],
                                      self.process_info['lon'],
                                      self.process_info['car_spd'])

        imps_trend = self.search_by_trend(self.process_info['time'],
                                          self.process_info['l_pos'],
                                          self.process_info['lat'],
                                          self.process_info['lon'],
                                          self.process_info['car_spd'])



        print("number of imps:", len(imps_accel), len(imps_spd), len(imps_spurt), len(imps_trend))




        fname = self.csvFileName.split('/')[-1]
        print("Filename is:", fname)
        fname_t = "accel_" + fname
        self.write_to_csv(fname_t, imps_accel)
        fname_t = "spd_" + fname
        self.write_to_csv(fname_t, imps_spd)
        fname_t = "spurt_" + fname
        self.write_to_csv(fname_t, imps_spurt)
        fname_t = "trend_" + fname
        self.write_to_csv(fname_t, imps_trend)


        filtered_data = self.filter_data(self.process_info['l_pos'],cutoff=2)
        plt.figure(1)
        ax1 = plt.subplot(411)
        plt.plot(self.process_info['time'], self.process_info['l_pos'], 'k')
        plt.plot(self.process_info['time'], filtered_data, 'g')
        plt.ylabel("Позиция, м")
        plt.grid(True)

        ax2 = plt.subplot(412, sharex=ax1)
        plt.plot(self.process_info['time'], self.process_info['l_spd'], 'k')
        plt.ylabel("Скорость, м/с")
        plt.grid(True)

        ax3 = plt.subplot(413, sharex=ax1)
        plt.plot(self.process_info['time'], self.process_info['l_accel'], 'k')
        plt.ylabel("Ускорение, м/с2")
        plt.grid(True)

        ax4 = plt.subplot(414, sharex=ax1)
        plt.plot(self.process_info['time'], self.process_info['l_spurt'], 'k')
        plt.ylabel("Рывок, м/с3")
        plt.grid(True)

        plt.show()

    def write_to_csv(self, csv_name, data):
        try:
            print(self.csvFileName)
            f = open(csv_name, 'w')
            f.write('Coordinates\n')

            for imp in data:
                if imp[1] != 0.0:
                    f.write('\"' + str(imp[1]) + ',' + str(imp[2]) + '\"\n')
            f.close()
        except Exception as ex:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)

    def filter_data(self, data, cutoff=0.5, fs=40, order=5):
        filtered_data = butter_lowpass_filter(data, cutoff, fs, order)

        tau = 1 / cutoff
        ztau = tau / 0.025 * 2

        for i in range(int(ztau)):
            filtered_data[i] = data[i]

        return filtered_data

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
    app = QApplication(sys.argv)
    ex = DataVisualisation()
    #sys.exit(app.exec_())