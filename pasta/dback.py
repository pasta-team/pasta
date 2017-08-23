# coding: utf-8
# Copyright © 2016 YunXing Zuo, WeiJi Hsiao

import numpy as np
import re
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import bisect

__author__ = 'JianShu Jie, WeiJi Hsiao, YunXing Zuo'
__email__ = 'weiji.hsiao@gmail.com'
__date__ = 'Mar. 16, 2017'

color = ['black','red','blue','lime','mediumvioletred','darkviolet','magenta','deepskyblue','springgreen','orangered','darkgoldenrod','orange','deepskyblue','peru','dodgerblue',\
'darkviolet','springgreen','darkgray','darkorange','chocolatesaddlebrown','dimgray','firebrick']


class Filedos(object):
    def __init__(self):
    	self.__energies = []
        self.__doses = []
        self.__efermi = None
        self.__titles = []

    def fermi_adjust(self):
        for i in range(len(self.__energies)):
            self.__energies[i] -= self.efermi

    def import_from_txt(self,dosfile,efermi = None):
        if len(self.__energies) == 0:
            spin = True
        else:
            spin = False
    	doslist = []
        if efermi is not None:
            self.__efermi = efermi
        flag = False
    	with open(dosfile) as f:
            while 1:
                line = f.readline()
                if not line:
                    break
                if not flag:
                    aline = re.sub('#','',line).strip()
                    self.__titles = aline.split()
                    flag = True
                else:
                    data = line.split()
                    for i,dosdata in enumerate(data):
                        dos = float(dosdata)
                        if i == 0:
                            if spin:
                                self.__energies.append(dos)
                        else:
                            doslist.append(dos)

        doses = np.array(doslist).reshape(-1,len(self.__titles)-1)
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


    def export_to_js(self,filename):
        with open(filename,'w') as f:
            for i in range(len(self.__doses[0])):
                f.write('var graphDos{:d} = '.format(i)+'{\n')
                f.write('data:[\n')
                for x,y in zip(self.energies,self.__doses[0][i]):
                    f.write('[{:8.4f},{:16.8f}],\n'.format(x,y))
                f.write(']};\n')
                if len(self.__doses) == 2:
                    f.write('var graphDos{:d} = '.format(i+self.no_doses)+'{\n')
                    f.write('data:[\n')
                    for x,y in zip(self.energies,self.__doses[1][i]):
                        f.write('[{:8.4f},{:16.8f}],\n'.format(x,y))
                    f.write(']};\n')
            f.write('var xmin = {:8.4f}, xmax = {:8.4f}, ymin = {:12.5e}, ymax = {:12.5e};\n'.format(self.emin,self.emax,self.dmin,self.dmax))
            f.write('var title = [')
            line = []
            for i in self.__titles:
                line.append("'"+i+"'")
            f.write(','.join(line))
            f.write('];')

    def export_for_analysis(self,filename):
        self.fermi_adjust()
        with open(filename,'w') as f:
            f.write('#')
            for i in self.titles:
                f.write('\t{:s}'.format(i))
            if len(self.__doses) == 2:
                for i in range(1,len(self.titles)):
                    f.write('\t{:s}'.format(self.titles[i]))
            f.write('\n')
            for i in range(len(self.energies)):
                f.write('{:8.4f}'.format(self.energies[i]))
                for j in range(len(self.__doses[0])):
                    f.write('\t{:12.5e}'.format(self.__doses[0][j][i]))
                if len(self.__doses) == 2:
                    for j in range(len(self.__doses[0])):
                        f.write('\t{:12.5e}'.format(self.__doses[1][j][i]))
                f.write('\n')

    '''
    def plot(self, foldname, funct, eps = True, svg = False):
        plt.figure()
        axes = plt.subplot(111)
        axes.set_xlim(self.energies[0],self.energies[-1])
        axes.set_ylim(1.2*self.dmin,1.2*self.dmax)
        axes.set_xlabel(r'$E-E_f\ (eV)$')

        base = len(self.titles)-1
        for i,y in enumerate(self.__doses[0]):
            plt.plot(self.energies,y,color[i%base],linewidth = 0.2,label = '{:s}'.format(self.titles[i+1]))
            plt.legend(loc = (1.01,1-base*0.028), fontsize = 5)
        if len(self.__doses) == 2:
            for i,y in enumerate(self.__doses[1]):
                plt.plot(self.energies,y,color[i%base],linewidth = 0.2,label = '{:s}'.format(self.titles[i+1]))
                plt.legend(loc = (1.01,1-2*base*0.028), fontsize = 5)

        if eps:
            plt.savefig('{:s}/{:s}dos.eps'.format(foldname, funct))
        if svg:
            plt.savefig('{:s}/{:s}dos.svg'.format(foldname, funct))
        plt.close()
    '''
    def plot(self, foldname, funct, eps = True, svg = False):
        plt.figure()
        axes = plt.subplot(111)
        emin = -5
        emax = 5
        l = bisect.bisect(self.energies,emin)-1
        r = bisect.bisect_left(self.energies,emax)+1
        x = self.energies[l:r+1]
        axes.set_xlim(self.energies[l],self.energies[r])
        ymin = float('inf')
        ymax = -float('inf')
        axes.set_ylim(1.2*self.dmin,1.2*self.dmax)
        axes.set_xlabel(r'$E-E_f\ (eV)$')

        base = len(self.titles)-1
        for i,y in enumerate(self.__doses[0]):
            d = y[l:r+1]
            ymax = max(ymax,max(d))
            plt.plot(x,d,color[i%base],linewidth = 0.2,label = '{:s}'.format(self.titles[i+1]))
            plt.legend(loc = (1.01,1-base*0.028), fontsize = 5)
        if len(self.__doses) == 2:
            for i,y in enumerate(self.__doses[1]):
                d = y[l:r+1]
                ymin = min(ymin,min(d))
                plt.plot(x,d,color[i%base],linewidth = 0.2,label = '{:s}'.format(self.titles[i+1]))
                plt.legend(loc = (1.01,1-2*base*0.028), fontsize = 5)
        else:
            ymin = 0
        axes.set_ylim(1.2*ymin,1.2*ymax)

        if eps:
            plt.savefig('{:s}/{:s}dos.eps'.format(foldname, funct))
        if svg:
            plt.savefig('{:s}/{:s}dos.svg'.format(foldname, funct))
        plt.close()

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
        if len(self.__doses) == 2:
            res = min(res,np.min(self.__doses[1]))
        return res

    @property
    def dmax(self):
        res = np.max(self.__doses[0])
        if len(self.__doses) == 2:
            res = max(res,np.max(self.__doses[1]))
        return res

    @property
    def no_doses(self):
        return self.__doses[0].shape[0]

    @property
    def dos_length(self):
        return self.__doses[0].shape[1]