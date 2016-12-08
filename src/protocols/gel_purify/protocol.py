from __future__ import print_function
from transcriptic_tools.utils import (ul, ml, ng,
                                      convert_mass_to_volume
                                      )

from transcriptic_tools.harness import run
from transcriptic_tools import CustomProtocol as Protocol
from transcriptic_tools.enums import Antibiotic, Reagent, Temperature

def get_purify(p, 
               dna_wells,
               measure_concentratio=True
            ):
    """
    
    Run a simple gel purification on one or more wells with dna
    
    """
    
    assert isinstance(p,Protocol)
   
    source_wells = []
    dna_lengths = []
    
    for dna_well in dna_wells:
        source_wells.append(dna_well)
        dna_lengths.append(int(dna_well.properties['dna_length']))
    
    dest_wells = p.simple_gel_purify(source_wells, dna_lengths)
   
    if measure_concentratio:
        p.measure_concentration(dest_wells, 'dna_concentration', 'DNA')
        
   

def main(p, params):    
    """This protocol takes a tube of bacteria, amplifies it, and creates/stores new tubes of frozen cells
    """
    #bacterial protocol
    p.mammalian_cell_mode = False
    
    get_purify(p,params['source_dna'],params['measure_concentration'])

if __name__ == '__main__':
    run(main, "GelPurify")
