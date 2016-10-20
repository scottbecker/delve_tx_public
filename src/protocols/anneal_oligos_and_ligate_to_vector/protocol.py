from __future__ import print_function
from transcriptic_tools.utils import (ul, ml, ng,
                                      convert_mass_to_volume
                                      )

from transcriptic_tools.harness import run
from transcriptic_tools import CustomProtocol as Protocol
from transcriptic_tools.enums import Antibiotic, Reagent, Temperature
import numpy

def anneal_oligos_and_ligate_to_vector(p, oligo1_tube, oligo2_tube,
                                       linearized_vector_tube,
                                       negative_control=True):
    """
    
    
    """
    
    
    #@TODO: update to dynamically detect concentration of primer tubes (currently assumes its 100uM)
    #@TODO: include option to proceed immediately to spread/transform/pick
    
    assert isinstance(p,Protocol)
    
    pcr_plate = p.ref('pcr_plate', cont_type="96-pcr", storage=Temperature.cold_4.name)
    anneal_well = pcr_plate.well(0)
    anneal_well.name = 'annealed_oligos'
    experiment_ligate_well = pcr_plate.well(1)
    experiment_ligate_well.name = 'ligated_vector_with_oligos'
    ligate_wells = [experiment_ligate_well]
    
    
    if negative_control:
        negative_control_well = pcr_plate.well(2)
        negative_control_well.name = 'negative_control_ligated_vector_only'
        ligate_wells.append(negative_control_well)
    
    atp_25mM_well = pcr_plate.well(3)
    atp_25mM_well.name = 'atp_25mM'
    
    p.provision_by_name(Reagent.water, atp_25mM_well, ul(6))
    p.provision_by_name(Reagent.atp_100mM, atp_25mM_well, ul(2))        
    
    #prep annealing, 10uL total
    
    p.provision_by_name(Reagent.water, anneal_well, ul(14))
    p.provision_by_name(Reagent.t4_polynucleotide_kinase_buffer_a_10x, anneal_well, ul(2))
    p.provision_by_name(Reagent.t4_polynucleotide_kinase, anneal_well, ul(1))
    #atp hasn't been mixed yet
    p.transfer(atp_25mM_well, anneal_well, ul(1), mix_before=True, mix_after=True)
    p.transfer([oligo1_tube.well(0),oligo2_tube.well(0)], anneal_well, ul(1),
               mix_after=True)
    
    
    #~45minutes = 38seconds * 70
    cooldown_steps = [{"temperature": "%s:celsius"%temp, 
                       "duration": 
                       "10:second"} for temp in numpy.arange(94,3.9,-1)]
    
    p.incubate(pcr_plate, Temperature.warm_37, '30:minute', shaking=False)
    
    p.thermocycle(pcr_plate, [{"cycles":  1, 
                               "steps": [{"temperature": "95:celsius", 
                                          "duration": 
                                          "5:minute"}
                                         ]+cooldown_steps}
                              ],volume=anneal_well.volume)
    
    #ligation in 20uL total volume
    
    vector_volume = convert_mass_to_volume(ng(20),linearized_vector_tube.well(0))
    
    p.provision_by_name(Reagent.water, ligate_wells, ul(16)-vector_volume)
    p.transfer(linearized_vector_tube.well(0), ligate_wells, vector_volume,
               mix_before=True, mix_after=True)    
    p.provision_by_name(Reagent.t4_ligase_kit_ligase_buffer_10x, ligate_wells, ul(2))
    p.transfer(anneal_well, experiment_ligate_well, ul(1), mix_after=True)    
    p.provision_by_name(Reagent.t4_ligase_kit_ligase, ligate_wells, ul(1),
                        mix_after=True)
    
    if negative_control:
        p.provision_by_name(Reagent.water, negative_control_well, ul(1), mix_after=True)

    
    p.incubate(pcr_plate, Temperature.ambient, '10:minute', shaking=False)
    
    #inactivate ligase
    p.thermocycle(pcr_plate, [{"cycles":  1, 
                               "steps": [
                                   {"temperature": "65:celsius", 
                                          "duration": 
                                          "10:minute"}
                                         ]}
                              ], volume=experiment_ligate_well.volume)
    
    

def main(p, params):    
    """This protocol takes a tube of bacteria, amplifies it, and creates/stores new tubes of frozen cells
    """
    #bacterial protocol
    p.mammalian_cell_mode = False
    
    anneal_oligos_and_ligate_to_vector(p,params['oligo1_tube'],params['oligo2_tube'],params['linearized_vector_tube'],
                                       negative_control=params['negative_control'])

if __name__ == '__main__':
    run(main, "AnnealOligosAndLigateToVector")
