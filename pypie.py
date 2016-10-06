import os
import numpy as np
import matplotlib.pyplot as plt
# import math

from collections import namedtuple
from uncertainties import unumpy

PieSpectrum = namedtuple('PieSpectrum', 'label center width counts')
MassSpectrum = namedtuple('MassSpectrum', 'label energy counts')

class Pypie(object):
    def __init__(self):
        """Returns a Pypie object with data from path.

        Args:
            path (str) : File path to PIE data.

        Returns:
            Pypie object: 
        """
        super().__init__()
        self.k, self.t0 = 1, 0
        self.pie_spectra = []
        self.mass_spectra = []
        self.data = {}
    def load(self, path):
        """Load PIE data from Oleg software.

        delimiter: \\t, newline: \\n, format:


        COL│    1      │     2     │   ∙∙∙   │    N    │ROW
        ───┴───────────┴───────────┴─────────┴─────────┼───
                                    "energy"           │ 
            ┌---------┐ ┌-----------------------------┐│
             'Energy'      float       ∙∙∙      float  │ 1
                                                       ├───
                                   "current"           │ 
            ┌---------┐ ┌-----------------------------┐│
             'Current'     float       ∙∙∙      float  │ 2
                                                       ├───
                               "data"                  │  
            ┌-----------------------------------------┐│ 
                int         int        ∙∙∙       int   │ 3
                                                       ├─── 
                 ∙           ∙          ∙         ∙    │ ∙
                 ∙           ∙          ∙         ∙    │ ∙
                 ∙           ∙          ∙         ∙    │ ∙
                                                       ├───
                int         int        ∙∙∙       int   │ N
            └---------┘ └-----------------------------┘│
               "time"              "counts"            │

        """
        self.path = path
        self.dir, self.name = os.path.split(path)
        data = np.genfromtxt(path)
        self.energy = np.around(data[0,1:], 3)
        self.current = data[1,1:]
        self.time = data[2:,0].astype(int)
        self.indices = np.asarray(range(self.time.size))
        self.mass = self.indices
        self.counts = data[2:,1:].astype(int)
        self.counts_sum = np.sum(self.counts, axis=1)
    def mass_calibrate(self, m1, m2, t1, t2):
        """ Calculate the proportionality, k, and time zero, t0, constants
        given two masses and their time/data points for t = k*sqrt(mass)

        Args:
            m1 (float) : mass 1 at time 1
            m2 (float) : mass 2 at time 2
            t1 (float) : time 1 at mass 1
            t2 (float) : time 2 at mass 2

        """
        self.t0 = -1 * ((-t2 + t1) \
            / (np.sqrt(m1) - np.sqrt(m2))) * np.sqrt(m1) + t1
        self.k = (-t2 + t1) / (np.sqrt(m1) - np.sqrt(m2))
        self.mass = [ self.index_to_mass(index) for index in self.indices ]
    def index_to_mass(self, index):
        """ Returns mass at a given point """
        return ((index-self.t0)/self.k)**2
    def mass_to_index(self, mass):
        """ Returns mass at a given point """
        return int(self.k*np.sqrt(mass) + self.t0)
    def pie_slice(self, center, width, label='', current=False):
        """ Sum over a subset of an """
        if not label: 
            label = ','.join([str(center), str(width)])
        half_width = width/2
        left, right = center - half_width, center + half_width
        indices = [ self.mass_to_index(i) for i in (left, right) ] + [1]
        indices = slice(*indices)
        counts = self.counts[indices, :]
        counts = np.sum(counts, axis=0)
        if current: 
            current = np.divide(self.current, np.amax(self.current))
            counts = np.divide(counts, current)
        pie_spectrum = PieSpectrum(label, center, width, counts)
        self.pie_spectra.append(pie_spectrum)
    def mass_slice(self, energy=None):
        """ Generate a mass obejct extracted from self.counts

        Args:
            energy (None) : Spectrum from integral over all counts  
            energy (float) : Spectrum at energy of counts
            energy (list[float,float]) : Spectrum from integral between energies 

        """
        if energy is None:
            energy = self.energy[0], self.energy[-1]
        elif type(energy) is type(0.):
            index = self.energy.searchsorted(energy)
            energy = energy, energy
        indices = slice(
            self.energy.searchsorted(energy[0]),
            self.energy.searchsorted(energy[1])+1,
            1)
        counts = self.counts[:,indices]
        counts = np.sum(counts, axis=1)
        label = '{}-{}'.format(*energy)
        mass_spectrum = MassSpectrum(label, energy, counts)
        self.mass_spectra.append(mass_spectrum)
    def energy_to_index(self, energy):
        """ Returns the index of energy from self.energy """
        return self.energy.searchsorted(energy)
    def mz_peak_edges(self, mass_spectrum, center, minimum):
        center = self.mass_to_index(center)
        left_min = center - np.argmax(mass_spectrum[:center][::-1] <= minimum) 
        right_min = center + np.argmax(mass_spectrum[center:] <= minimum)
        width = self.mass[right_min] - self.mass[left_min]
        return width

class uPypie(Pypie):
    """Returns a Pypie object with propgating uncertainty from data in paths.

    paths : list, tuple
        File paths to n PIE data sets

    """
    def __init__(self):
        super().__init__()
    def load(self, paths):
        """Load multiple PIE data sets with uncertainty using Oleg format 

        See Pypie.load for data format.

        """
        self.paths = paths
        counts,  energy,  current = [], [], []
        for path in paths:
            data = np.genfromtxt(path)
            energy.append(np.around(data[0,1:], 3))
            current.append(data[1,1:])
            counts.append(data[2:,1:].astype(int))
        avg, std = np.average(energy, axis=0), np.std(energy, axis=0)
        self.energy = unumpy.uarray(avg, std)
        avg, std = np.average(current, axis=0), np.std(current, axis=0)
        self.current = unumpy.uarray(avg, std)
        avg, std = np.average(counts, axis=0), np.std(counts, axis=0)
        self.counts = unumpy.uarray(avg, std)
        # self.counts_sum = np.sum(self.counts, axis=1)
        self.time = data[2:,0].astype(int)
        self.indices = np.asarray(range(self.time.size))
        self.mass = self.indices
    def pie_save(self, path):
        # header = 'Energy\t\t' 'Counts'*len(data)
        header = [ spectrum.label for spectrum in self.pie_spectra ]
        header = '\t±\t'.join(['eV'] + header) + '\t±'
        data = [ spectrum.counts for spectrum in self.pie_spectra ]
        fmt = ('%.3f',)*2
        fmt += ('%d',)*2 * len(data)
        # data = [ [ unumpy.nominal_values(ar), unumpy.std_devs(ar) ] \
            # for ar in data ]
        # data = np.vstack(data).astype(int)
        # energy, uncertainty = unumpy.nominal_values(self.energy), unumpy.std_devs(self.energy)
        data = np.vstack([ self.split_val_stdv(ar) for ar in data ])
        data = data.astype(int)
        energy, uncertainty = self.split_val_stdv(self.energy)
        data = np.vstack([energy, uncertainty, data]).T
        np.savetxt(
            path, 
            data, 
            fmt=fmt, 
            header=header, 
            comments='',
            delimiter='\t')
        # np.savetxt(path, data, fmt='%r', header=header, comments='')
    def current_save(self, path):
        header = 'Energy /eV\tCurrent /A\n'
        data = np.vstack([ self.split_val_stdv(ar) \
            for ar in [ self.energy,  self.current] ])
        data = None
    @staticmethod
    def split_val_stdv(uarray):
        """ Returns the nominal values and uncertainties as separate arrays. """
        # values, uncertainties = np.array([ (x.n, x.s) for x in uarray ]).T
        values = unumpy.nominal_values(uarray)
        uncertainties = unumpy.std_devs(uarray)
        return values, uncertainties

def pie_plot(pypie_obj):
    ax = plt.gca()
    energy = unumpy.nominal_values(pypie_obj.energy)
    pie_spectra = pypie_obj.pie_spectra
    for spectrum in pie_spectra:
        counts = unumpy.nominal_values(spectrum.counts)
        ax.plot(energy, counts, label=spectrum.label)
    plt.legend()
def mass_plot(pypie_obj, energy):
    ax = plt.gca()
    mass = pypie_obj.mass
    mass_spectra = pypie_obj.pie_spectra

if __name__ == '__main__':
    pass