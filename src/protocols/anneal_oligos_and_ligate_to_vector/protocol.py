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
                                       negative_control=True,
                                       transform_spread_pick=False,
                                       antibiotic=None):
    """
    
    
    """
    
    
    #@TODO: update to dynamically detect concentration of primer tubes (currently assumes its 100uM)
    #@TODO: include option to proceed immediately to spread/transform/pick
    
    assert isinstance(p,Protocol)
    
    reaction_plate = p.ref('reaction_plate', cont_type="96-pcr", storage=Temperature.cold_4.name)
    anneal_well = reaction_plate.well(0)
    anneal_well.name = 'annealed_oligos'

    
    atp_25mM_well = reaction_plate.well('C1')
    atp_25mM_well.name = 'atp_25mM'
    
    p.provision_by_name(Reagent.water, atp_25mM_well, ul(6))
    p.provision_by_name(Reagent.atp_100mM, atp_25mM_well, ul(2))        
    
    oligo_wells = [oligo1_tube.well(0),oligo2_tube.well(0)]
    
    assert all([well.properties.get('Concentration')=='100uM' for well in oligo_wells]), 'All oligos must be at 100uM'
    
    #------------- Anneal / Phosphorylation ----------------------------
    
    #prep annealing, 20uL total
    #see instructions on p1 here: https://www.addgene.org/static/data/plasmids/60/60958/60958-attachment_wWVpb-8u9Mzp.pdf
    
    p.provision_by_name(Reagent.water, anneal_well, ul(14))
    p.provision_by_name(Reagent.t4_polynucleotide_kinase_buffer_a_10x, anneal_well, ul(2))
    p.provision_by_name(Reagent.t4_polynucleotide_kinase, anneal_well, ul(1))
    #atp hasn't been mixed yet
    p.transfer(atp_25mM_well, anneal_well, ul(1), mix_before=True, mix_after=True)
    p.transfer(oligo_wells, anneal_well, ul(1),
               mix_after=True)
    
    
    assert anneal_well.volume == ul(20)
    
    p.incubate(reaction_plate, Temperature.warm_37, '30:minute', shaking=False)
    
    #~45minutes = 38seconds * 70
    cooldown_steps = [{"temperature": "%s:celsius"%temp, 
                       "duration": 
                       "10:second"} for temp in numpy.arange(94,3.9,-1)]    
    
    p.thermocycle(reaction_plate, [{"cycles":  1, 
                               "steps": [{"temperature": "95:celsius", 
                                          "duration": 
                                          "5:minute"}
                                         ]+cooldown_steps}
                              ],volume=anneal_well.volume)
    
    # ------------ Dilute ----------------
    
    dilution_wells = reaction_plate.wells_from('H1', 3)
    
    #10x,100x,500x
    
    
    dilution_wells[0].name = '10x_oligo_dilution'    
    dilution_wells[1].name = '100x_oligo_dilution'    
    dilution_wells[2].name = '500x_oligo_dilution'

    p.provision_by_name(Reagent.water, dilution_wells[:2], ul(45))
    p.provision_by_name(Reagent.water, dilution_wells[2],ul(20))
    
    previous_well = anneal_well
    for x in range(0,3):
        dilution_well = dilution_wells[x]
        p.transfer(previous_well, dilution_well, ul(5), mix_after=True)
        
        previous_well = dilution_well
    
        
    #------------- Ligate ----------------------------
    

    experiment_ligate_wells = reaction_plate.wells_from(1,len(dilution_wells))
    for i, experiment_ligate_well in enumerate(experiment_ligate_wells):
        experiment_ligate_well.name = 'ligated_vector_with_%s'%dilution_wells[i].name
        
    ligate_wells = list(experiment_ligate_wells)
    
    if negative_control:
        negative_control_well = reaction_plate.well('B1')
        negative_control_well.name = 'negative_control_ligated_vector_only'
        ligate_wells.append(negative_control_well)    
    
    
    #ligation in 20uL total volume
    total_volume = ul(20)
    vector_volume = convert_mass_to_volume(ng(20),linearized_vector_tube.well(0))
    t4_buffer_volume = ul(2)
    t4_ligase_volume = ul(1)
    oligo_volume = ul(1)
    
    water_volume = total_volume-vector_volume-t4_buffer_volume-t4_ligase_volume-oligo_volume
    
    p.provision_by_name(Reagent.water, ligate_wells,water_volume )
    p.transfer(linearized_vector_tube.well(0), ligate_wells, vector_volume,
               mix_before=True, mix_after=True)    
    p.provision_by_name(Reagent.t4_ligase_kit_ligase_buffer_10x, ligate_wells, t4_buffer_volume)
    
    p.transfer(dilution_wells, experiment_ligate_wells, oligo_volume, mix_after=True)    
    
    if negative_control:
        p.provision_by_name(Reagent.water, negative_control_well, oligo_volume, mix_after=True)
    
    
    #ligase added last
    p.provision_by_name(Reagent.t4_ligase_kit_ligase, ligate_wells, t4_ligase_volume,
                        mix_after=True)
    
    
    assert all([well.volume == total_volume for well in ligate_wells]), 'all ligation wells must be 20uL'
    
    
    p.incubate(reaction_plate, Temperature.ambient, '10:minute', shaking=False)
    
    #inactivate ligase
    p.thermocycle(reaction_plate, [{"cycles":  1, 
                               "steps": [
                                   {"temperature": "65:celsius", 
                                          "duration": 
                                          "10:minute"}
                                         ]}
                              ], volume=ligate_wells[0].volume)
    
    #transform, spread, pick 2 (recommended to follow up with sequencing)
    #we have 4 wells to spread
    if transform_spread_pick:
        p.transform_spread_pick(ligate_wells, antibiotic,
                                last_source_well_is_negative_control=True)
    
    

def main(p, params):    
    """This protocol takes a tube of bacteria, amplifies it, and creates/stores new tubes of frozen cells
    """
    #bacterial protocol
    p.mammalian_cell_mode = False
    
    anneal_oligos_and_ligate_to_vector(p,params['oligo1_tube'],params['oligo2_tube'],params['linearized_vector_tube'],
                                       negative_control=params['negative_control'],
                                       transform_spread_pick=params['transform_spread_pick'],
                                       antibiotic=Antibiotic.from_string(params['antibiotic']),
                                       )

if __name__ == '__main__':
    run(main, "AnnealOligosAndLigateToVector")
