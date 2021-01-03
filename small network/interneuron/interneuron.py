from neuron import h
from neuron.units import ms,mV
from cell import Cell
class Interneuron(Cell):
    name = 'Interneuron'
    def _set_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.soma.L=self.soma.diam = 10
        
    def _set_biophysics(self):

        for sec in self.all:
            sec.Ra = 100
            sec.cm = 1.5
            
        

        