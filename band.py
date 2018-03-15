# coding: utf-8
# Copyright © 2016 YunXing Zuo, WeiJi Hsiao

import numpy as np
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

__author__ = 'JianShu Jie, WeiJi Hsiao, YunXing Zuo'
__email__ = 'weiji.hsiao@gmail.com'
__date__ = 'Mar. 16, 2017'


class Fileband(object):
    def __init__(self):
        pass
        self.__x = []
        self.__path = []
        self.__high_symmetry_points = []
        self.__bands = []
        self.__efermi = None

    def fermi_adjust(self):
        for i in self.__bands:
            i -= self.__efermi
        self.__efermi = 0

    def is_metal(self):
        for i in range(self.no_bands):
            if np.any(self.__bands[0][i, :] < self.__efermi) and np.any(self.__bands[0][i, :] > self.efermi):
                return True
            if len(self.__bands) == 2:
                if np.any(self.__bands[1][i, :] < self.__efermi) and np.any(self.__bands[1][i, :] > self.efermi):
                    return True
        return False

    def vbm(self):
        if self.is_metal():
            return None
        else:
            res = max(list(v for v in self.__bands[0].flat if v < self.efermi))
            if len(self.__bands) == 2:
                res = max(
                    res, max(list(v for v in self.__bands[1].flat if v < self.efermi)))
            return res

    def cbm(self):
        if self.is_metal():
            return None
        else:
            res = min(list(v for v in self.__bands[0].flat if v > self.efermi))
            if len(self.__bands) == 2:
                res = min(
                    res, min(list(v for v in self.__bands[1].flat if v > self.efermi)))
            return res

    def band_gap(self):
        if self.is_metal():
            return 0.0
        else:
            return self.cbm() - self.vbm()

    def import_from_txt(self, bandfile, path, efermi=None):
        if len(self.__x) == 0:
            spin = True
        else:
            spin = False
        if spin:
            self.__path = path
        if efermi is not None:
            self.__efermi = efermi
        flag = False
        bandlist = []
        no_bands = 0
        with open(bandfile) as f:
            while 1:
                line = f.readline()
                if not line:
                    break
                if re.match('\s+$', line):
                    no_bands += 1
                    if flag == False:
                        flag = True
                else:
                    x, band_point = map(float, line.split())
                    bandlist.append(band_point)
                    if spin and not flag:
                        self.__x.append(x)
        self.__bands.append(np.array(bandlist).reshape(no_bands, -1))
        if spin:
            for (no, label) in path:
                self.__high_symmetry_points.append((self.__x[int(no)], label))

    def export_to_js(self, filename):
        with open(filename, 'w') as f:
            f.write('var graphBand = [\n')
            for band in self.__bands[0]:
                f.write('{data:[\n')
                for x, y in zip(self.x, band):
                    f.write('[{:8.4f},{:16.8f}],\n'.format(x, y))
                f.write("],\ncolor:'#000000',},\n")
            if len(self.__bands) == 2:
                for band in self.__bands[1]:
                    f.write('{data:[\n')
                    for x, y in zip(self.x, band):
                        f.write('[{:8.4f},{:16.8f}],\n'.format(x, y))
                    f.write("],\ncolor:'#000000',},\n")
            f.write('];\nvar xmin = {:8.4f}, xmax = {:8.4f}, ymin = {:16.10f}, ymax = {:16.10f};'.format(
                self.xmin, self.xmax, self.bmin, self.bmax))
            f.write('\nvar HighSymmetryPoints = [')
            for x, y in self.high_symmetry_points:
                f.write(" [{:8.4f}, '{:s}'],".format(x, y))
            f.write('];')
            if len(self.__bands) == 2:
                f.write("\nvar spin = 'True'")
            else:
                f.write("\nvar spin = 'False'")

    def export_for_analysis(self, filename):
        if len(self.__bands) == 1:
            f = open(filename, 'w')
            for band in self.__bands[0]:
                for x, y in zip(self.x, band):
                    f.write('{:8.4f}\t{:16.8f}\n'.format(x, y))
                f.write('\n')
        else:
            f = open('{:s}spin'.format(filename), 'w')
            for band1, band2 in zip(self.__bands[0], self.__bands[1]):
                for x, y1, y2 in zip(self.x, band1, band2):
                    f.write('{:8.4f}\t{:16.8f}\t{:16.8f}\n'.format(x, y1, y2))
                f.write('\n')
        f.close()

    def plot(self, foldname, funct, eps=True, png=False):
        plt.figure()
        axes = plt.subplot(111)
        axes.set_xlim(0, self.x[-1])
        emin = self.efermi - 5
        emax = self.efermi + 5
        axes.set_ylim(emin, emax)
        axes.set_xticks([])
        axes.set_yticks([emin, (emin + emax) / 2, emax])
        axes.set_ylabel(r'$E-E_f\ (eV)$')
        axes.set_xlabel('Wave Vector')

        for i in self.__bands[0]:
            for point in i:
                if point < emin or point > emax:
                    break
            else:
                plt.plot(self.x, i, '-b', linewidth=0.2)

        if len(self.__bands) == 2:
            for i in self.__bands[1]:
                for point in i:
                    if point < emin or point > emax:
                        break
                else:
                    plt.plot(self.x, i, '-r', linewidth=0.2)

        z = []
        for i in self.x:
            z.append(self.efermi)

        plt.plot(self.x, z, '--k', linewidth=0.2)

        xt1 = []
        xt2 = []
        for (i, j) in self.high_symmetry_points:
            plt.axvline(float(i), ls='-', color='k', linewidth=0.2)
            xt1.append(float(i))
            xt2.append(r'${:s}$'.format(j))
            #plt.annotate(r'${:s}$'.format(j),xy = (float(i),-5.5),color = 'k',size = 5)

        tag = 0
        for i in range(len(xt1) - 1):
            if xt1[i] == xt1[i + 1]:
                xt1.append(xt1[i])
                xt1[i] -= 0.03
                xt1[i + 1] += 0.03
                tag += 1
                xt2.insert(i + tag, r'$|$')
        xt1.sort()
        axes.tick_params(axis='x', length=0)
        plt.xticks(xt1, xt2, size=8)
        plt.yticks(size=8)
        if eps:
            plt.savefig('{:s}/{:s}band.eps'.format(foldname, funct))
        if png:
            plt.savefig('{:s}/{:s}band.png'.format(foldname, funct), dpi=600)
        plt.close()

    def plot_with_fermi(self, foldname, funct, eps=True, png=False):
        plt.figure()
        axes = plt.subplot(111)
        axes.set_xlim(0, self.x[-1])
        emin = self.efermi - 5
        emax = self.efermi + 5
        axes.set_ylim(emin, emax)
        axes.set_xticks([])
        axes.set_yticks([emin, (emin + emax) / 2, emax])
        axes.set_ylabel(r'$E\ (eV)$')
        axes.set_xlabel('Wave Vector')

        for i in self.__bands[0]:
            for point in i:
                if point < emin or point > emax:
                    break
            else:
                plt.plot(self.x, i, '-b', linewidth=0.2)

        if len(self.__bands) == 2:
            for i in self.__bands[1]:
                for point in i:
                    if point < emin or point > emax:
                        break
                else:
                    plt.plot(self.x, i, '-r', linewidth=0.2)

        z = []
        for i in self.x:
            z.append(self.efermi)

        plt.plot(self.x, z, '--k', linewidth=0.2)

        xt1 = []
        xt2 = []
        for (i, j) in self.high_symmetry_points:
            plt.axvline(float(i), ls='-', color='k', linewidth=0.2)
            xt1.append(float(i))
            xt2.append(r'${:s}$'.format(j))
            #plt.annotate(r'${:s}$'.format(j),xy = (float(i),-5.5),color = 'k',size = 5)

        tag = 0
        for i in range(len(xt1) - 1):
            if xt1[i] == xt1[i + 1]:
                xt1.append(xt1[i])
                xt1[i] -= 0.03
                xt1[i + 1] += 0.03
                tag += 1
                xt2.insert(i + tag, r'$|$')
        xt1.sort()
        axes.tick_params(axis='x', length=0)
        plt.xticks(xt1, xt2, size=8)
        plt.yticks(size=8)
        if eps:
            plt.savefig('{:s}/{:s}band.eps'.format(foldname, funct))
        if png:
            plt.savefig('{:s}/{:s}band.png'.format(foldname, funct), dpi=600)
        plt.close()

    @property
    def efermi(self):
        return self.__efermi

    @property
    def x(self):
        return self.__x

    @property
    def path(self):
        return self.__path

    @property
    def high_symmetry_points(self):
        return self.__high_symmetry_points

    @property
    def no_bands(self):
        return self.__no_bands

    @property
    def band_length(self):
        return self.__band_length

    @property
    def xmin(self):
        return np.min(self.x)

    @property
    def xmax(self):
        return np.max(self.x)

    @property
    def bmin(self):
        return np.min(self.__bands)

    @property
    def bmax(self):
        return np.max(self.__bands)

    @property
    def no_bands(self):
        return self.__bands[0].shape[0]

    @property
    def band_length(self):
        return self.__bands[0].shape[1]
