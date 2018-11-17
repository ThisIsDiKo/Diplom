from PyQt5.QtWidgets import QFileDialog, QMainWindow, QApplication
import sys
import matplotlib.pyplot as plt
from statistics import median, variance, mean
from math import sqrt
import numpy as np
from scipy.signal import butter, lfilter, freqz


class CalibrationModel(QMainWindow):
    def __init__(self):
        super().__init__()

        self.timings_up = []
        self.timings_down = []

        self.l_pos_up = []
        self.l_press_up = []
        self.r_pos_up = []
        self.r_press_up = []

        self.l_pos_down = []
        self.l_press_down= []
        self.r_pos_down = []
        self.r_press_down = []

        self.suspension_model = {'max': 0,
                                 'min': 0,
                                 'border_up': 0,
                                 'border_down': 0,
                                 'k1': 0,
                                 'r1': 0,
                                 'a1': 0,
                                 'b1': 0,
                                 'c1': 0,
                                 'a2': 0,
                                 'b2': 0,
                                 'c2': 0,
                                 'a3': 0,
                                 'b3': 0,
                                 'c3': 0}

        self.open_file()
        #self.invert_results()
        self.simulate_calibration_one_circuit()
    def invert_results(self):
        r = 2480
        for i in range(len(self.l_pos_up)):
            self.l_pos_up[i] = r + (r - self.l_pos_up[i])
        for i in range(len(self.l_pos_down)):
            self.l_pos_down[i] = r + (r - self.l_pos_down[i])

    def open_file(self):
        fname_down = QFileDialog.getOpenFileName(self, 'Движение вниз', 'C:\\Users\\ADiKo\\Desktop\\demos', 'txt file (*.txt)')[0]
        fname_up = QFileDialog.getOpenFileName(self, 'Движение вверх', 'C:\\Users\\ADiKo\\Desktop\\demos', 'txt file (*.txt)')[0]

        try:
            f = open(fname_down, 'r')
            for line in f:
                if line and len(line) > 5:
                    l = line.split(',')
                    if l[0] is 'a':
                        self.timings_down.append(int(l[1]))
                        self.l_pos_down.append(int(l[2]))
                        self.l_press_down.append(int(l[3]))
                        self.r_pos_down.append(int(l[4]))
                        self.r_press_down.append(int(l[5]))
            f.close()

            f = open(fname_up, 'r')
            for line in f:
                if line and len(line) > 5:
                    l = line.split(',')
                    if l[0] is 'a':
                        self.timings_up.append(int(l[1]))
                        self.l_pos_up.append(int(l[2]))
                        self.l_press_up.append(int(l[3]))
                        self.r_pos_up.append(int(l[4]))
                        self.r_press_up.append(int(l[5]))
            f.close()
            t = self.timings_down[0]
            for i in range(len(self.timings_down)):
                self.timings_down[i] -= t
            t = self.timings_up[0]
            for i in range(len(self.timings_up)):
                self.timings_up[i] -= t

        except Exception as ex:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)

    def simulate_calibration_one_circuit(self):
        NUM_OF_INTERVALS = 64
        threshold = int(3000 / NUM_OF_INTERVALS)
        #threshold = 50
        mass = []
        prev_pos = self.l_pos_down[0]
        for p in self.l_pos_down:
            if abs(prev_pos - p) > threshold:
                mass.append([self.timings_down[self.l_pos_down.index(p)],
                             p,
                             self.l_press_down[self.l_pos_down.index(p)]])
                prev_pos = p



        temp = int(len(mass) / 3)

        (k1, r1) = get_linear_coeff(mass[0][0], mass[0][1],
                                    mass[len(mass) - temp][0],  mass[len(mass) - temp][1])
        (a1, b1, c1) = get_parabola_coeff(mass[len(mass) - temp][0],  mass[len(mass) - temp][1],
                                          mass[int(len(mass) - temp/2)][0], mass[int(len(mass) - temp/2)][1],
                                          mass[len(mass) - 1][0], mass[len(mass) - 1][1])

        prev_pos = self.l_pos_up[0]
        mass_up = []
        for p in self.l_pos_up:
            if abs(prev_pos - p) > threshold:
                mass_up.append([self.timings_up[self.l_pos_up.index(p)],
                                p,
                                self.l_press_up[self.l_pos_up.index(p)]])
                prev_pos = p



        temp_up= int(len(mass_up) / 2)  # вообще 3, но пох)
        (a2, b2, c2) = get_parabola_coeff(mass_up[0][0], mass_up[0][1],
                                          mass_up[int(temp_up / 2)][0], mass_up[int(temp_up / 2)][1],
                                          mass_up[temp_up][0], mass_up[temp_up][1])
        (a3, b3, c3) = get_parabola_coeff(mass_up[temp_up][0], mass_up[temp_up][1],
                                          mass_up[int(len(mass_up) - temp_up / 2)][0], mass_up[int(len(mass_up) - temp_up / 2)][1],
                                          mass_up[len(mass_up) - 1][0], mass_up[len(mass_up) - 1][1])


        self.suspension_model['max'] = self.l_pos_down[0]
        self.suspension_model['min'] = self.l_pos_up[0]
        self.suspension_model['border_down'] = mass[len(mass) - temp][1]
        self.suspension_model['border_up'] = mass_up[temp_up][1]
        self.suspension_model['k1'] = k1
        self.suspension_model['r1'] = r1
        self.suspension_model['a1'] = a1
        self.suspension_model['b1'] = b1
        self.suspension_model['c1'] = c1
        self.suspension_model['a2'] = a2
        self.suspension_model['b2'] = b2
        self.suspension_model['c2'] = c2
        self.suspension_model['a3'] = a3
        self.suspension_model['b3'] = b3
        self.suspension_model['c3'] = c3
        if self.suspension_model['max'] > self.suspension_model['min']:
            self.suspension_model['inverted'] = False
        else:
            self.suspension_model['inverted'] = True


        print(self.suspension_model)

        dot_x = []
        dot_y = []
        dot_x.append(self.get_t_up(2200))
        dot_y.append(2200)
        dot_x.append(self.get_t_up(2900))
        dot_y.append(2900)
        dot_x.append(self.get_t_up(2600))
        dot_y.append(2600)

        dot_x_d = []
        dot_y_d = []
        dot_x_d.append(self.get_t_down(2200))
        dot_y_d.append(2200)
        dot_x_d.append(self.get_t_down(2900))
        dot_y_d.append(2900)
        dot_x_d.append(self.get_t_down(2600))
        dot_y_d.append(2600)




        down_approx = {'x': [], 'y': []}
        up_approx = {'x': [], 'y': []}

        q = self.l_pos_down.index(mass[0][1])
        w = self.l_pos_down.index(mass[len(mass) - temp][1])
        x1 = self.timings_down[q:w+1]
        y1 = [k1 * i + r1 for i in x1]

        q = w + 1
        w = self.l_pos_down.index(mass[len(mass) - 1][1])
        x2 = self.timings_down[q:w+1]
        y2 = [a1 * i * i + b1 * i + c1 for i in x2]
        down_approx['x'].extend(x1)
        down_approx['x'].extend(x2)
        down_approx['y'].extend(y1)
        down_approx['y'].extend(y2)

        q = self.l_pos_up.index(mass_up[0][1])
        w = self.l_pos_up.index(mass_up[temp_up][1])
        x1 = self.timings_up[q:w + 1]
        y1 = [a2 * i * i + b2 * i + c2 for i in x1]

        q = self.l_pos_up.index(mass_up[temp_up][1]) + 1
        w = self.l_pos_up.index(mass_up[len(mass_up) - 1][1])
        x2 = self.timings_up[q:w + 1]
        y2 = [a3 * i * i + b3 * i + c3 for i in x2]

        up_approx['x'].extend(x1)
        up_approx['x'].extend(x2)
        up_approx['y'].extend(y1)
        up_approx['y'].extend(y2)

        plt.figure(1)
        plt.plot(self.timings_down, self.l_pos_down, 'r')
        plt.plot(down_approx['x'], down_approx['y'], 'b')
        plt.scatter(dot_x_d, dot_y_d, color="black")
        plt.figure(2)
        plt.plot(self.timings_up, self.l_pos_up, 'r')
        plt.plot(up_approx['x'], up_approx['y'], 'b')
        #plt.plot(dot_x, dot_y,  marker='o', markersize=6, color="black")
        plt.scatter(dot_x, dot_y,color="black")

        plt.grid(True)
        plt.show()

    def get_t_up(self, pos):
        root = 0
        if pos < self.suspension_model['border_up']:
            if self.suspension_model['inverted']:
                D = self.suspension_model['b3'] * self.suspension_model['b3'] - 4 * self.suspension_model['a3'] * (
                            self.suspension_model['c3'] - pos)
                root1 = (-self.suspension_model['b3'] + sqrt(D)) / (2 * self.suspension_model['a3'])
                root2 = (-self.suspension_model['b3'] - sqrt(D)) / (2 * self.suspension_model['a3'])
                print("i roots:", root1, root2)
                if root1 > root2:
                    return root2
                else:
                    return root1
            else:
                D = self.suspension_model['b2'] * self.suspension_model['b2'] - 4 * self.suspension_model['a2'] * (
                            self.suspension_model['c2'] - pos)

                root1 = (-self.suspension_model['b2'] + sqrt(D)) / (2 * self.suspension_model['a2'])
                root2 = (-self.suspension_model['b2'] - sqrt(D)) / (2 * self.suspension_model['a2'])
                print("roots:", root1, root2)
                if root1 > root2:
                    return root1
                else:
                    return root2
        else:
            if self.suspension_model['inverted']:
                D = self.suspension_model['b2'] * self.suspension_model['b2'] - 4 * self.suspension_model['a2'] * (
                            self.suspension_model['c2'] - pos)
                root1 = (-self.suspension_model['b2'] + sqrt(D)) / (2 * self.suspension_model['a2'])
                root2 = (-self.suspension_model['b2'] - sqrt(D)) / (2 * self.suspension_model['a2'])
                print("i roots:", root1, root2)
                if root1 > root2:
                    return root1
                else:
                    return root2
            else:
                D = self.suspension_model['b3'] * self.suspension_model['b3'] - 4 * self.suspension_model['a3'] * (
                            self.suspension_model['c3'] - pos)
                root1 = (-self.suspension_model['b3'] + sqrt(D)) / (2 * self.suspension_model['a3'])
                root2 = (-self.suspension_model['b3'] - sqrt(D)) / (2 * self.suspension_model['a3'])
                print("roots:", root1, root2)
                if root1 > root2:
                    return root2
                else:
                    return root1


    def get_t_down(self, pos):
        if pos < self.suspension_model['border_down']:
            if self.suspension_model['inverted']:
                root = (pos - self.suspension_model['r1']) / self.suspension_model['k1']
                print("i root:", root)
                return root
            else:
                D = self.suspension_model['b1'] * self.suspension_model['b1'] - 4 * self.suspension_model['a1'] * (
                        self.suspension_model['c1'] - pos)

                root1 = (-self.suspension_model['b1'] + sqrt(D)) / (2 * self.suspension_model['a1'])
                root2 = (-self.suspension_model['b1'] - sqrt(D)) / (2 * self.suspension_model['a1'])
                print("roots:", root1, root2)
                if root1 < root2:
                    return root1
                else:
                    return root2
        else:
            if self.suspension_model['inverted']:
                D = self.suspension_model['b1'] * self.suspension_model['b1'] - 4 * self.suspension_model['a1'] * (
                        self.suspension_model['c1'] - pos)
                root1 = (-self.suspension_model['b1'] + sqrt(D)) / (2 * self.suspension_model['a1'])
                root2 = (-self.suspension_model['b1'] - sqrt(D)) / (2 * self.suspension_model['a1'])
                print("i roots:", root1, root2)
                if root1 < root2:
                    return root1
                else:
                    return root2
            else:
                root = (pos - self.suspension_model['r1']) / self.suspension_model['k1']
                print("root:", root)
                return root

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


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = CalibrationModel()