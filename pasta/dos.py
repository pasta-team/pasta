# coding: utf-8
# Copyright © 2016 YunXing Zuo, WeiJi Hsiao

import numpy as np
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import bisect

__author__ = 'JianShu Jie, WeiJi Hsiao, YunXing Zuo'
__email__ = 'weiji.hsiao@gmail.com'
__date__ = 'Mar. 16, 2017'

# color = ['black', 'red', 'blue', 'lime', 'mediumvioletred', 'darkviolet', 'magenta', 'deepskyblue', 'springgreen', 'orangered', 'darkgoldenrod', 'orange', 'deepskyblue', 'peru', 'dodgerblue', 'darkviolet', 'springgreen', 'darkgray', 'darkorange', 'chocolatesaddlebrown', 'dimgray', 'firebrick']

color = ['#000000','#FF0000','#9AFF9A','#FFFF00','#00FFFF','#9400D3','#FFA500','#EEAEEE','#D2691E','#006400','#0000FF','#F5F5DC','#F4A460','#C0C0C0','#FA8072','#FFD700']


class Filedos(object):
    def __init__(self):
        self.__energies = []
        self.__doses = []
        self.__efermi = None
        self.__titles = []
        self.__inte = {}
        self.__efermi_backup = None

    def fermi_adjust(self):
        for i in range(len(self.__energies)):
            self.__energies[i] -= self.efermi
        self.__efermi_backup = self.efermi
        self.__efermi = 0
    
    def fermi_restore(self):
        self.__efermi = self.__efermi_backup
        for i in range(len(self.__energies)):
            self.__energies[i] += self.efermi

    def import_from_txt(self, dosfile, efermi=None):
        if len(self.__energies) == 0:
            spin = True
        else:
            spin = False
        doslist = []
        if efermi is not None:
            self.__efermi = efermi
        flag = False
        anoflag = False
        with open(dosfile) as f:
            while 1:
                line = f.readline()
                if not line:
                    break
                if not flag:
                    aline = re.sub('#', '', line).strip()
                    self.__titles = aline.split()
                    flag = True
                else:
                    data = line.split()
                    if not anoflag:
                        anoflag = True
                        if len(data) != len(self.__titles):
                            raise ValueError
                    for i, dosdata in enumerate(data):
                        dos = float(dosdata)
                        if i == 0:
                            if spin:
                                self.__energies.append(dos)
                        else:
                            if not spin and dos > 0:
                                doslist.append(-dos)
                            else:
                                doslist.append(dos)

        doses = np.array(doslist).reshape(-1, len(self.__titles) - 1)
        doses = np.transpose(doses)
        self.__doses.append(doses)

    '''
    def export_to_csv(self,filename):
    	with open(filename,'w') as f:
            f.write('#')
            for i in self.titles:
                f.write('\t{:s}'.format(i))
            f.write('\n')
            for i in range(len(self.energies)):
                f.write('\"{:8.4f}\"'.format(self.energies[i]))
                for j in range(len(self.__doses)):
                    f.write(',\"{:12.5e}\"'.format(self.__doses[j][i]))
                f.write('\n')
    '''

    def export_to_js(self, filename):
        with open(filename, 'w') as f:
            for i in range(len(self.__doses[0])):
                f.write('var graphDos{:d} = '.format(i) + '{\n')
                f.write('data:[\n')
                for x, y in zip(self.energies, self.__doses[0][i]):
                    f.write('[{:8.4f},{:16.8f}],\n'.format(x, y))
                f.write(']};\n')
                if self.spin:
                    f.write('var graphDos{:d} = '.format(
                        i + self.no_doses) + '{\n')
                    f.write('data:[\n')
                    for x, y in zip(self.energies, self.__doses[1][i]):
                        f.write('[{:8.4f},{:16.8f}],\n'.format(x, y))
                    f.write(']};\n')
            f.write('var xmin = {:8.4f}, xmax = {:8.4f}, ymin = {:12.5e}, ymax = {:12.5e};\n'.format(
                self.emin, self.emax, self.dmin, self.dmax))
            f.write('var title = [')
            line = []
            for i in self.__titles:
                line.append("'" + i + "'")
            f.write(','.join(line))
            f.write('];')

    def export_for_analysis(self, filename):
        with open(filename, 'w') as f:
            f.write('#')
            for i in self.titles:
                f.write('\t{:s}'.format(i))
            if self.spin:
                for i in range(1, len(self.titles)):
                    f.write('\t{:s}'.format(self.titles[i]))
            f.write('\n')
            for i in range(len(self.energies)):
                f.write('{:8.4f}'.format(self.energies[i]))
                for j in range(len(self.__doses[0])):
                    f.write('\t{:12.5e}'.format(self.__doses[0][j][i]))
                if self.spin:
                    for j in range(len(self.__doses[0])):
                        f.write('\t{:12.5e}'.format(self.__doses[1][j][i]))
                f.write('\n')

    def plot(self, foldname, funct, eps=True, png=False, erela=True, emin=None, emax=None, elimit = 0.03):
        plt.figure()
        axes = plt.subplot(111)
        if erela:
            if emin is None:
                emin = -5 + self.efermi
            else:
                emin = emin + self.efermi
            if emax is None:
                emax = 5 + self.efermi
            else:
                emax = emax + self.efermi
        else:
            if emin is None or emax is None:
                raise ValueError('Energy range should be specified in absolute mode')
        l = bisect.bisect(self.energies, emin)
        if l > 0:
            l -= 1
        r = bisect.bisect_left(self.energies, emax)
        if r < len(self.energies) - 1:
            r += 1
        if r >= len(self.energies):
            r -= 1
        x = self.energies[l:r + 1]
        axes.set_xlim(self.energies[l], self.energies[r])
        ymin = float('inf')
        ymax = -float('inf')
        axes.set_ylim(1.2 * self.dmin, 1.2 * self.dmax)
        axes.set_xlabel(r'$E\ (eV)$')
        plt.axvline(self.efermi, ls='--', color='k', linewidth=0.2)
        z = []
        for i in x:
            z.append(0)

        plt.plot(x, z, '--k', linewidth=0.2)

        for dos in self.__doses[0]:
            d = dos[l:r+1]
            ymax = max(ymax, max(d))
        elimit1 = ymax * elimit
        if self.spin:
            for dos in self.__doses[1]:
                d = dos[l:r+1]
                ymin = min(ymin, min(d))
            elimit2 = ymin * elimit
        else:
            ymin = 0

        i = 0
        ploted = {}
        base = len(self.titles) - 1
        for j, y in enumerate(self.__doses[0]):
            d = y[l:r + 1]
            if max(d) > elimit1:
                plt.plot(x, d, color[i], linewidth=0.2,
                         label='{:s}'.format(self.titles[j + 1]))
                plt.legend(loc=(1.01, 1 - base * 0.028), fontsize=5)
                ploted[j] = i
                i += 1
        if self.spin:
            i = len(ploted)
            for j, y in enumerate(self.__doses[1]):
                d = y[l:r + 1]
                if min(d) < elimit2:
                    if j in ploted:
                        plt.plot(x, d, color[ploted[j]], linewidth=0.2,label='{:s}'.format(self.titles[j + 1]))
                    else:
                        plt.plot(x, d, color[i], linewidth=0.2,label='{:s}'.format(self.titles[j + 1]))
                        plt.legend(loc=(1.01, 1 - base * 0.028), fontsize=5)
                        i += 1

        axes.set_ylim(1.2 * ymin, 1.2 * ymax)

        if eps:
            plt.savefig('{:s}/{:s}dos.eps'.format(foldname, funct))
        if png:
            plt.savefig('{:s}/{:s}dos.png'.format(foldname, funct), dpi=600)
        plt.close()

    def plot_backup(self, foldname, funct, eps=True, png=False):
        plt.figure()
        axes = plt.subplot(111)
        emin = -5
        emax = 5
        l = bisect.bisect(self.energies, emin)
        if l > 0:
            l -= 1
        r = bisect.bisect_left(self.energies, emax)
        if r < len(self.energies) - 1:
            r += 1
        if r >= len(self.energies):
            r -= 1
        x = self.energies[l:r + 1]
        axes.set_xlim(self.energies[l], self.energies[r])
        ymin = float('inf')
        ymax = -float('inf')
        axes.set_xlabel(r'$E-E_f\ (eV)$')

        base = len(self.titles) - 1
        for i, y in enumerate(self.__doses[0]):
            d = y[l:r + 1]
            ymax = max(ymax, max(d))
            plt.plot(x, d, color[i % base], linewidth=0.2,
                     label='{:s}'.format(self.titles[i + 1]))
            plt.legend(loc=(1.01, 1 - base * 0.028), fontsize=5)
        if self.spin:
            for i, y in enumerate(self.__doses[1]):
                d = y[l:r + 1]
                ymin = min(ymin, min(d))
                plt.plot(x, d, color[i % base], linewidth=0.2,
                         label='{:s}'.format(self.titles[i + 1]))
                plt.legend(loc=(1.01, 1 - 2 * base * 0.028), fontsize=5)
        else:
            ymin = 0
        if ymin == ymax:
            print(foldname)
        axes.set_ylim(1.2 * ymin, 1.2 * ymax)

        if eps:
            plt.savefig('{:s}/{:s}dos.eps'.format(foldname, funct))
        if png:
            plt.savefig('{:s}/{:s}dos.png'.format(foldname, funct), dpi=600)
        plt.close()
    
    def dos_inte(self):
        l12 = bisect.bisect_right(self.energies,self.efermi-12)
        l9 = bisect.bisect_right(self.energies,self.efermi-9)
        l6 = bisect.bisect_right(self.energies,self.efermi-6)
        l3 = bisect.bisect_right(self.energies,self.efermi-3)
        l = bisect.bisect_right(self.energies,self.efermi)
        self.__inte['s']=[0]*4
        self.__inte['p']=[0]*4
        self.__inte['d']=[0]*4
        self.__inte['sum']=[0]*4
        sset = set()
        pset = set()
        dset = set()
        for i in range(2,len(self.titles)):
            data = self.titles[i].split('-')
            if data[1] == 's':
                sset.add(i)
            elif data[1] == 'p':
                pset.add(i)
            elif data[1] == 'd':
                dset.add(i)
        
        for i in range(l12,l9):
            self.__inte['sum'][0] += self.__doses[0][0][i]
            for j in range(2,len(self.titles)):
                if j in sset:
                    self.__inte['s'][0] += self.__doses[0][j-1][i]
                elif j in pset:
                    self.__inte['p'][0] += self.__doses[0][j-1][i]
                elif j in dset:
                    self.__inte['d'][0] += self.__doses[0][j-1][i]
        if self.spin:
            for i in range(l12,l9):
                self.__inte['sum'][0] -= self.__doses[1][0][i]
                for j in range(2,len(self.titles)):
                    if j in sset:
                        self.__inte['s'][0] -= self.__doses[1][j-1][i]
                    elif j in pset:
                        self.__inte['p'][0] -= self.__doses[1][j-1][i]
                    elif j in dset:
                        self.__inte['d'][0] -= self.__doses[1][j-1][i]
        
        for i in range(l9,l6):
            self.__inte['sum'][1] += self.__doses[0][0][i]
            for j in range(2,len(self.titles)):
                if j in sset:
                    self.__inte['s'][1] += self.__doses[0][j-1][i]
                elif j in pset:
                    self.__inte['p'][1] += self.__doses[0][j-1][i]
                elif j in dset:
                    self.__inte['d'][1] += self.__doses[0][j-1][i]
        if self.spin:
            for i in range(l9,l6):
                self.__inte['sum'][1] -= self.__doses[1][0][i]
                for j in range(2,len(self.titles)):
                    if j in sset:
                        self.__inte['s'][1] -= self.__doses[1][j-1][i]
                    elif j in pset:
                        self.__inte['p'][1] -= self.__doses[1][j-1][i]
                    elif j in dset:
                        self.__inte['d'][1] -= self.__doses[1][j-1][i]
        
        for i in range(l6,l3):
            self.__inte['sum'][2] += self.__doses[0][0][i]
            for j in range(2,len(self.titles)):
                if j in sset:
                    self.__inte['s'][2] += self.__doses[0][j-1][i]
                elif j in pset:
                    self.__inte['p'][2] += self.__doses[0][j-1][i]
                elif j in dset:
                    self.__inte['d'][2] += self.__doses[0][j-1][i]
        if self.spin:
            for i in range(l6,l3):
                self.__inte['sum'][2] -= self.__doses[1][0][i]
                for j in range(2,len(self.titles)):
                    if j in sset:
                        self.__inte['s'][2] -= self.__doses[1][j-1][i]
                    elif j in pset:
                        self.__inte['p'][2] -= self.__doses[1][j-1][i]
                    elif j in dset:
                        self.__inte['d'][2] -= self.__doses[1][j-1][i]
        
        for i in range(l3,l):
            self.__inte['sum'][3] += self.__doses[0][0][i]
            for j in range(2,len(self.titles)):
                if j in sset:
                    self.__inte['s'][3] += self.__doses[0][j-1][i]
                elif j in pset:
                    self.__inte['p'][3] += self.__doses[0][j-1][i]
                elif j in dset:
                    self.__inte['d'][3] += self.__doses[0][j-1][i]
        if self.spin:
            for i in range(l3,l):
                self.__inte['sum'][3] -= self.__doses[1][0][i]
                for j in range(2,len(self.titles)):
                    if j in sset:
                        self.__inte['s'][3] -= self.__doses[1][j-1][i]
                    elif j in pset:
                        self.__inte['p'][3] -= self.__doses[1][j-1][i]
                    elif j in dset:
                        self.__inte['d'][3] -= self.__doses[1][j-1][i]
    
    @property
    def spd_percent(self):
        li = []
        for i in range(4):
            s = self.__inte['s'][i]+self.__inte['p'][i]+self.__inte['d'][i]
            if s < 1e-5:
                li.append([0,0,0])
            else:
                li.append([self.__inte['s'][i]/s,self.__inte['p'][i]/s,self.__inte['d'][i]/s])
        return li
    
    @property
    def dos_percent(self):
        s = [0]*4
        for i in range(4):
            s[i] = self.__inte['s'][i]+self.__inte['p'][i]+self.__inte['d'][i]
        ss = sum(s)
        li = []
        for i in range(4):
            li.append(s[i]/ss)
        return li


    @property
    def energies(self):
        return self.__energies

    @property
    def titles(self):
        return self.__titles

    @property
    def efermi(self):
        return self.__efermi

    @property
    def emin(self):
        return np.min(self.energies)

    @property
    def emax(self):
        return np.max(self.energies)

    @property
    def dmin(self):
        res = np.min(self.__doses[0])
        if self.spin:
            res = min(res, np.min(self.__doses[1]))
        return res

    @property
    def dmax(self):
        res = np.max(self.__doses[0])
        if self.spin:
            res = max(res, np.max(self.__doses[1]))
        return res

    @property
    def no_doses(self):
        return self.__doses[0].shape[0]

    @property
    def dos_length(self):
        return self.__doses[0].shape[1]
    
    @property
    def dos_at_fermi(self):
        l = bisect.bisect_left(self.energies, self.efermi)
        if self.efermi-self.energies[l] > self.energies[l+1]-self.efermi:
            l += 1
        daf = self.__doses[0][0][l]
        if self.spin:
            daf -= self.__doses[1][0][l]
        return daf
    
    @property
    def spin(self):
        return len(self.__doses) == 2
