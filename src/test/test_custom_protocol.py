from __future__ import print_function
import unittest
from datadiff import diff
from transcriptic_tools.utils import (ul, ml, get_well_max_volume, get_well_dead_volume,
                                      hours, get_column_wells)
from test_utils import create_blank_plate
from transcriptic_tools.custom_protocol import CustomProtocol as Protocol
from autoprotocol.instruction import Provision, Incubate, Cover, Dispense, Absorbance
from transcriptic_tools.enums import Reagent, Antibiotic

class TestCustomProtocol(unittest.TestCase):
    maxDiff = None
    
    #@TODO add version that uses shared medium or freezing
    
    def test_culture_medium_provisioning(self):
        
        initial_protocol = Protocol(mammalian_cell_mode=True)
        
        plate = initial_protocol.ref("my_plate", cont_type="24-deep", 
                                     storage="cold_4", discard=False)
        
        initial_protocol.provision_by_name(Reagent.culture_medium, plate.well(0), ml(1),
                                           mix_after=True)
        
        p = initial_protocol._get_final_protocol()
        
        #check that the instructions include provisioning statements at the beginning for
        assert isinstance(p.instructions[0],Provision)
        
        self.assertEqual(p.instructions[0].data['resource_id'],
                         p.transcriptic_inv[Reagent.dmem_fbs_ps])

        #check that there is a cover statement before the first incubate of culture medium
        
        
        for i in range(0,len(p.instructions)):
            instruction = p.instructions[i]
            if isinstance(instruction,Incubate):
                self.assertEqual(type(p.instructions[i-1]),Cover)
        
        #final protocol needs to be the same as the initial protocol
        self.assertEqual(p.as_dict(seal_on_store=False),initial_protocol.as_dict(seal_on_store=False))
        self.assertDictEqual(p.as_dict(seal_on_store=False),
                             {
                                 'refs': {
                                     'my_plate': {
                                         'new': '24-deep',
                                         'store': {
                                             'where': 'cold_4'
                                         }
                                     },
                                     'culture_medium': {
                                         'new': '96-deep',
                                         'discard': True
                                     }
                                 },
                                 'time_constraints': [{
                                     'to': {
                                         'instruction_start': 3
                                     },
                                     'less_than': '1.0:hour',
                                     'from': {
                                         'instruction_end': 1
                                     }
                                 }, {
                                     'to': {
                                         'instruction_end': 1
                                     },
                                     'less_than': '1.0:hour',
                                     'from': {
                                         'instruction_start': 3
                                     }
                                 }],
                                 'instructions': [{
                                     'to': [{
                                         'volume': '900.0:microliter',
                                         'well': 'culture_medium/0'
                                     }, {
                                         'volume': '115.0:microliter',
                                         'well': 'culture_medium/0'
                                     }],
                                     'op': 'provision',
                                     'resource_id': 'rs197gzgq2fufr'
                                 }, {
                                     'lid': 'standard',
                                     'object': 'culture_medium',
                                     'op': 'cover'
                                 }, {
                                     'where': 'warm_37',
                                     'object': 'culture_medium',
                                     'co2_percent': 0,
                                     'duration': '5:minute',
                                     'shaking': False,
                                     'op': 'incubate'
                                 }, {
                                     'object': 'culture_medium',
                                     'op': 'uncover'
                                 }, {
                                     'groups': [{
                                         'transfer': [{
                                             'volume': '900.0:microliter',
                                             'to': 'my_plate/0',
                                             'from': 'culture_medium/0',
                                             'mix_before': {
                                                 'volume': '900.0:microliter',
                                                 'repetitions': 5,
                                                 'speed': '900.0:microliter/second'
                                             }
                                         }, {
                                             'volume': '100.0:microliter',
                                             'to': 'my_plate/0',
                                             'from': 'culture_medium/0',
                                             'mix_after': {
                                                 'volume': '500.0:microliter',
                                                 'repetitions': 5,
                                                 'speed': '500.0:microliter/second'
                                             }
                                         }]
                                     }],
                                     'op': 'pipette'
                                 }]
                             }                         
                         )
        
    def test_culture_medium_provisioning_filled_dest(self):
            
        initial_protocol = Protocol(mammalian_cell_mode=True)
            
        plate = initial_protocol.ref("my_plate", cont_type="24-deep", 
                                     storage="cold_4", discard=False)
        
        for well in plate.all_wells():
            well.volume = ul(100)
        
        initial_protocol.provision_by_name('culture_medium', plate.well(0), ml(1))
        
        p = initial_protocol._get_final_protocol()
        
        #check that the instructions include provisioning statements at the beginning for
        assert isinstance(p.instructions[0],Provision)
        
        
        #check that the instructions include provisioning statements at the beginning for
        assert isinstance(p.instructions[0],Provision)
   
        self.assertEqual(p.instructions[0].data['resource_id'],
                        p.transcriptic_inv[Reagent.dmem_fbs_ps])
   
        #check that there is a cover statement before the first incubate of culture medium
       
   
        for i in range(0,len(p.instructions)):
            instruction = p.instructions[i]
            if isinstance(instruction,Incubate):
                self.assertEqual(type(p.instructions[i-1]),Cover)
            
        #final protocol needs to be the same as the initial protocol
        self.assertEqual(p.as_dict(seal_on_store=False),initial_protocol.as_dict(seal_on_store=False))
        
        expected_result =  {
            'refs': {
                'my_plate': {
                    'new': '24-deep',
                    'store': {
                        'where': 'cold_4'
                    }
                },
                'culture_medium': {
                    'new': '96-deep',
                    'discard': True
                }
            },
            'time_constraints': [{
                'to': {
                    'instruction_start': 3
                },
                'less_than': '1.0:hour',
                'from': {
                    'instruction_end': 1
                }
            }, {
                'to': {
                    'instruction_end': 1
                },
                'less_than': '1.0:hour',
                'from': {
                    'instruction_start': 3
                }
            }],
            'instructions': [{
                'to': [{
                    'volume': '900.0:microliter',
                    'well': 'culture_medium/0'
                }, {
                    'volume': '115.0:microliter',
                    'well': 'culture_medium/0'
                }],
                'op': 'provision',
                'resource_id': 'rs197gzgq2fufr'
            }, {
                'lid': 'standard',
                'object': 'culture_medium',
                'op': 'cover'
            }, {
                'where': 'warm_37',
                'object': 'culture_medium',
                'co2_percent': 0,
                'duration': '5:minute',
                'shaking': False,
                'op': 'incubate'
            }, {
                'object': 'culture_medium',
                'op': 'uncover'
            }, {
                'groups': [{
                    'transfer': [{
                        'volume': '900.0:microliter',
                        'to': 'my_plate/0',
                        'from': 'culture_medium/0',
                        'mix_before': {
                            'volume': '900.0:microliter',
                            'repetitions': 5,
                            'speed': '900.0:microliter/second'
                        }
                    }]
                }, {
                    'transfer': [{
                        'volume': '100.0:microliter',
                        'to': 'my_plate/0',
                        'from': 'culture_medium/0'
                    }]
                }],
                'op': 'pipette'
            }]
        }
        
        self.assertDictEqual(p.as_dict(seal_on_store=False),expected_result
                                                    
                    )    
    def test_culture_medium_large_provisioning(self):
    
        initial_protocol = Protocol()
    
        plate = initial_protocol.ref("my_plate", cont_type="24-deep", 
                                     storage="cold_4", discard=False)
    
        initial_protocol.provision_by_name(Reagent.culture_medium, plate.all_wells(), ml(5),
                                           mix_after=True)
    
        p = initial_protocol._get_final_protocol()
    
        #check that the instructions include provisioning statements at the beginning for
        assert isinstance(p.instructions[0],Dispense)
        
        
        self.assertEqual(p.instructions[0].data['resource_id'],
                         p.transcriptic_inv[Reagent.dmem_fbs_ps])
    
        #check that there is a cover statement before the first incubate of culture medium
    
        for i in range(0,len(p.instructions)):
            instruction = p.instructions[i]
            if isinstance(instruction,Incubate):
                self.assertEqual(type(p.instructions[i-1]),Cover)
    
        protocol_dict = p.as_dict(seal_on_store=False)    
    
        #final protocol needs to be the same as the initial protocol
        self.assertEqual(protocol_dict,initial_protocol.as_dict(seal_on_store=False))
        self.assertDictEqual({
            'object': 'culture_medium',
            'op': 'dispense',
            'columns': [{
                'column': 0,
                'volume': '2000.0:microliter'
            }, {
                'column': 1,
                'volume': '2000.0:microliter'
            }, {
                'column': 2,
                'volume': '2000.0:microliter'
            }, {
                'column': 3,
                'volume': '2000.0:microliter'
            }, {
                'column': 4,
                'volume': '2000.0:microliter'
            }, {
                'column': 5,
                'volume': '2000.0:microliter'
            }, {
                'column': 6,
                'volume': '2000.0:microliter'
            }, {
                'column': 7,
                'volume': '1120.0:microliter'
            }],
            'resource_id': 'rs197gzgq2fufr'
        },
        protocol_dict['instructions'][0])
        
    def test_column_provision_converted_to_dispense(self):
        
        initial_protocol = Protocol()
    
        plate = initial_protocol.ref("my_plate", cont_type="96-deep", 
                                     storage="cold_4", discard=False)
    
        initial_protocol.provision_by_name(Reagent.pbs, get_column_wells(plate,[0,11]), ml(1),
                                           mix_after=True)
    
        p = initial_protocol._get_final_protocol()
    
        #check that the instructions include provisioning statements at the beginning for
        assert isinstance(p.instructions[0],Dispense)
        
        self.assertEqual(p.instructions[0].data['resource_id'],
                         p.transcriptic_inv[Reagent.pbs])
    
        protocol_dict = p.as_dict(seal_on_store=False)    
    
        #final protocol needs to be the same as the initial protocol
        self.assertDictEqual({
            'object': 'my_plate',
            'op': 'dispense',
            'columns': [{
                'column': 0,
                'volume': '1000.0:microliter'
            }, {
                'column': 11,
                'volume': '1000.0:microliter'
            }],
            'resource_id': 'rs194na2u3hfam'
        },
        protocol_dict['instructions'][0])
        
       

    def test_remove_redundant_mixing_before(self):
        p = Protocol()
        
        transfer_instruction_group = [{
            'transfer': [{
                    'volume': '900.0:microliter',
                    'to': 'my_plate/0',
                    'from': 'culture_medium/0',
                    'mix_before': {
                        'volume': '900.0:microliter',
                        'repetitions': 5,
                        'speed': '900.0:microliter/second'
                    }
                }, {
                    'volume': '100.0:microliter',
                    'to': 'my_plate/0',
                    'from': 'culture_medium/0',
                    'mix_before': {
                        'volume': '900.0:microliter',
                        'repetitions': 5,
                        'speed': '900.0:microliter/second'
                    }
                }]
            }]
        
        p._remove_redundant_mixing(transfer_instruction_group)
        
        self.assertEqual(transfer_instruction_group,
                         [{
                                     'transfer': [{
                                         'volume': '900.0:microliter',
                                         'to': 'my_plate/0',
                                         'from': 'culture_medium/0',
                                         'mix_before': {
                                             'volume': '900.0:microliter',
                                             'repetitions': 5,
                                             'speed': '900.0:microliter/second'
                                         }
                                         }, {
                                             'volume': '100.0:microliter',
                                             'to': 'my_plate/0',
                                             'from': 'culture_medium/0'
                                         }]
                         }]
                         )
        
    def test_remove_redundant_mixing_before_include_distribute(self):
        p = Protocol()
        
        transfer_instruction_group = [{
            'transfer': [{
                    'volume': '900.0:microliter',
                    'to': 'my_plate/0',
                    'from': 'culture_medium/0',
                    'mix_before': {
                        'volume': '900.0:microliter',
                        'repetitions': 5,
                        'speed': '900.0:microliter/second'
                    }
                }, {
                    'volume': '100.0:microliter',
                    'to': 'my_plate/0',
                    'from': 'culture_medium/0',
                    'mix_before': {
                        'volume': '900.0:microliter',
                        'repetitions': 5,
                        'speed': '900.0:microliter/second'
                    }
                }]
            },
            {'distribute':{}},
            {
            'transfer': [{
                    'volume': '900.0:microliter',
                    'to': 'my_plate/0',
                    'from': 'culture_medium/0',
                    'mix_before': {
                        'volume': '900.0:microliter',
                        'repetitions': 5,
                        'speed': '900.0:microliter/second'
                    }
                }, {
                    'volume': '100.0:microliter',
                    'to': 'my_plate/0',
                    'from': 'culture_medium/0',
                    'mix_before': {
                        'volume': '900.0:microliter',
                        'repetitions': 5,
                        'speed': '900.0:microliter/second'
                    }
                }]
            }            
                                      
            ]
        
        p._remove_redundant_mixing(transfer_instruction_group)
        
        self.assertEqual(transfer_instruction_group,
                         [{
                                     'transfer': [{
                                         'volume': '900.0:microliter',
                                         'to': 'my_plate/0',
                                         'from': 'culture_medium/0',
                                         'mix_before': {
                                             'volume': '900.0:microliter',
                                             'repetitions': 5,
                                             'speed': '900.0:microliter/second'
                                         }
                                         }, {
                                             'volume': '100.0:microliter',
                                             'to': 'my_plate/0',
                                             'from': 'culture_medium/0'
                                         }]
                         },
                          {'distribute':{}},
                          
                          {
                              'transfer': [{
                                  'volume': '900.0:microliter',
                                  'to': 'my_plate/0',
                                  'from': 'culture_medium/0',
                                  'mix_before': {
                                      'volume': '900.0:microliter',
                                      'repetitions': 5,
                                      'speed': '900.0:microliter/second'
                                  }
                                  }, {
                                      'volume': '100.0:microliter',
                                      'to': 'my_plate/0',
                                      'from': 'culture_medium/0'
                                  }]
                              },                          
                          
                          
                          ]
                         )
    def test_remove_redundant_mixing_before_diff_dests(self):
        p = Protocol()
        
        transfer_instruction_group = [{
            'transfer': [{
                    'volume': '900.0:microliter',
                    'to': 'my_plate/0',
                    'from': 'culture_medium/0',
                    'mix_before': {
                        'volume': '900.0:microliter',
                        'repetitions': 5,
                        'speed': '900.0:microliter/second'
                    }
                }, {
                    'volume': '100.0:microliter',
                    'to': 'my_plate/0',
                    'from': 'culture_medium/0',
                    'mix_before': {
                        'volume': '900.0:microliter',
                        'repetitions': 5,
                        'speed': '900.0:microliter/second'
                    }
                }]
            },
              {
                  'transfer': [{
                      'volume': '900.0:microliter',
                      'to': 'my_plate/1',
                      'from': 'culture_medium/0',
                      'mix_before': {
                          'volume': '900.0:microliter',
                          'repetitions': 5,
                          'speed': '900.0:microliter/second'
                      }
                  }, {
                      'volume': '100.0:microliter',
                      'to': 'my_plate/1',
                      'from': 'culture_medium/0',
                      'mix_before': {
                          'volume': '900.0:microliter',
                          'repetitions': 5,
                          'speed': '900.0:microliter/second'
                      }
                  }]
              }
            ]     
        
        p._remove_redundant_mixing(transfer_instruction_group)
        
        self.assertEqual(transfer_instruction_group,
                         [{
                             'transfer': [{
                                 'volume': '900.0:microliter',
                                 'to': 'my_plate/0',
                                 'from': 'culture_medium/0',
                                 'mix_before': {
                                     'volume': '900.0:microliter',
                                     'repetitions': 5,
                                     'speed': '900.0:microliter/second'
                                 }
                                 }, {
                                     'volume': '100.0:microliter',
                                     'to': 'my_plate/0',
                                     'from': 'culture_medium/0'
                                 }]
                         },
                          
                          {
                              'transfer': [{
                                  'volume': '900.0:microliter',
                                  'to': 'my_plate/1',
                                  'from': 'culture_medium/0',                            
                              }, {
                                  'volume': '100.0:microliter',
                                  'to': 'my_plate/1',
                                  'from': 'culture_medium/0',
                              }]
                          }                          
                        ]

                         )    
           
    def test_remove_redundant_mixing_before_and_after(self):
        p = Protocol()
        
        transfer_instruction_group = [{
            'transfer': [{
                    'volume': '900.0:microliter',
                    'to': 'my_plate/0',
                    'from': 'culture_medium/0',
                    'mix_before': {
                        'volume': '900.0:microliter',
                        'repetitions': 5,
                        'speed': '900.0:microliter/second'
                    },
                    'mix_after': {
                        'volume': '500.0:microliter',
                        'repetitions': 5,
                        'speed': '500.0:microliter/second'
                    }                        
                }, {
                    'volume': '100.0:microliter',
                    'to': 'my_plate/0',
                    'from': 'culture_medium/0',
                    'mix_before': {
                        'volume': '100.0:microliter',
                        'repetitions': 5,
                        'speed': '100.0:microliter/second'
                    },
                    'mix_after': {
                        'volume': '500.0:microliter',
                        'repetitions': 5,
                        'speed': '500.0:microliter/second'
                    }                        
                }]
            }]
        
        p._remove_redundant_mixing(transfer_instruction_group)
        
        self.assertEqual(transfer_instruction_group,
                         [{
                                     'transfer': [{
                                         'volume': '900.0:microliter',
                                         'to': 'my_plate/0',
                                         'from': 'culture_medium/0',
                                         'mix_before': {
                                             'volume': '900.0:microliter',
                                             'repetitions': 5,
                                             'speed': '900.0:microliter/second'
                                         }
                                         }, {
                                             'volume': '100.0:microliter',
                                             'to': 'my_plate/0',
                                             'from': 'culture_medium/0',
                                             'mix_after': {
                                                 'volume': '500.0:microliter',
                                                 'repetitions': 5,
                                                 'speed': '500.0:microliter/second'
                                             }                                                   
                                         }]
                         }]
                         )    
        
        
    def test_remove_redundant_mixing_before_and_after_simple_mode(self):
            p = Protocol()
            
            transfer_instruction_group = [{
                'transfer': [{
                        'volume': '900.0:microliter',
                        'to': 'my_plate/0',
                        'from': 'culture_medium/0',
                        'mix_before': {
                            'volume': '900.0:microliter',
                            'repetitions': 5,
                            'speed': '900.0:microliter/second'
                        },
                        'mix_after': {
                            'volume': '500.0:microliter',
                            'repetitions': 5,
                            'speed': '500.0:microliter/second'
                        }                        
                    },
                            {
                                'volume': '900.0:microliter',
                                'to': 'my_plate/1',
                                'from': 'culture_medium/0',
                                'mix_before': {
                                    'volume': '900.0:microliter',
                                    'repetitions': 5,
                                    'speed': '900.0:microliter/second'
                                },
                                'mix_after': {
                                    'volume': '900.0:microliter',
                                    'repetitions': 5,
                                    'speed': '900.0:microliter/second'
                                }                                
                            }                             
                        , {
                        'volume': '100.0:microliter',
                        'to': 'my_plate/0',
                        'from': 'culture_medium/0',
                        'mix_before': {
                            'volume': '100.0:microliter',
                            'repetitions': 5,
                            'speed': '100.0:microliter/second'
                        },
                        'mix_after': {
                            'volume': '500.0:microliter',
                            'repetitions': 5,
                            'speed': '500.0:microliter/second'
                        }
                        },
                            {
                                'volume': '900.0:microliter',
                                'to': 'my_plate/1',
                                'from': 'culture_medium/1',
                                'mix_after': {
                                    'volume': '900.0:microliter',
                                    'repetitions': 5,
                                    'speed': '900.0:microliter/second'
                                }                     
                            }                       
                    ]
                }]
            
            p._remove_redundant_mixing(transfer_instruction_group)
            
            self.assertEqual(transfer_instruction_group,
                             [{
                                         'transfer': [{
                                             'volume': '900.0:microliter',
                                             'to': 'my_plate/0',
                                             'from': 'culture_medium/0',
                                             'mix_before': {
                                                 'volume': '900.0:microliter',
                                                 'repetitions': 5,
                                                 'speed': '900.0:microliter/second'
                                             }
                                             },{
                                                 'volume': '900.0:microliter',
                                                 'to': 'my_plate/1',
                                                 'from': 'culture_medium/0'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'to': 'my_plate/0',
                                                 'from': 'culture_medium/0',
                                                 'mix_after': {
                                                     'volume': '500.0:microliter',
                                                     'repetitions': 5,
                                                     'speed': '500.0:microliter/second'
                                                 }                                                   
                                             },
                                                    
                                             {
                                                 'volume': '900.0:microliter',
                                                 'to': 'my_plate/1',
                                                 'from': 'culture_medium/1',
                                                 'mix_after': {
                                                     'volume': '900.0:microliter',
                                                     'repetitions': 5,
                                                     'speed': '900.0:microliter/second'
                                                 }                     
                                             }                                                      
                             ]}]
                             )        

    def test_remove_redundant_mixing_before_multi_tip(self):
        p = Protocol()
        
        transfer_instruction_groups = [{
            'transfer': [{
                'volume': '900.0:microliter',
                'to': 'my_plate/0',
                'from': 'culture_medium/0',
                'mix_before': {
                    'volume': '900.0:microliter',
                    'repetitions': 5,
                    'speed': '900.0:microliter/second'
                }
            }]
        }, {
            'transfer': [{
                'volume': '100.0:microliter',
                'to': 'my_plate/0',
                'from': 'culture_medium/0',
                'mix_before': {
                    'volume': '900.0:microliter',
                    'repetitions': 5,
                    'speed': '900.0:microliter/second'
                }
            }]
        }]
        
        p._remove_redundant_mixing(transfer_instruction_groups)
        
        self.assertEqual(transfer_instruction_groups,
                         [{
                             'transfer': [{
                                 'volume': '900.0:microliter',
                                 'to': 'my_plate/0',
                                 'from': 'culture_medium/0',
                                 'mix_before': {
                                     'volume': '900.0:microliter',
                                     'repetitions': 5,
                                     'speed': '900.0:microliter/second'
                                 }
                             }]
                         }, {
                             'transfer': [{
                                 'volume': '100.0:microliter',
                                 'to': 'my_plate/0',
                                 'from': 'culture_medium/0'
                             }]
                         }]
                         )
        
    def test_remove_redundant_mixing_before_multi_tip_simple_mode(self):
        p = Protocol()
        
        transfer_instruction_groups = [{
            'transfer': [{
                'volume': '900.0:microliter',
                'to': 'my_plate/0',
                'from': 'culture_medium/0',
                'mix_before': {
                    'volume': '900.0:microliter',
                    'repetitions': 5,
                    'speed': '900.0:microliter/second'
                }
            },
            {
                'volume': '900.0:microliter',
                'to': 'my_plate/1',
                'from': 'culture_medium/0',
                'mix_before': {
                    'volume': '900.0:microliter',
                    'repetitions': 5,
                    'speed': '900.0:microliter/second'
                },
                'mix_after': {
                    'volume': '900.0:microliter',
                    'repetitions': 5,
                    'speed': '900.0:microliter/second'
                },                            
                
            }                         
            ]
        }, {
            'transfer': [{
                'volume': '100.0:microliter',
                'to': 'my_plate/0',
                'from': 'culture_medium/0',
                'mix_before': {
                    'volume': '900.0:microliter',
                    'repetitions': 5,
                    'speed': '900.0:microliter/second'
                }
            },
                         
            {
                'volume': '900.0:microliter',
                'to': 'my_plate/1',
                'from': 'culture_medium/1',
                'mix_after': {
                    'volume': '900.0:microliter',
                    'repetitions': 5,
                    'speed': '900.0:microliter/second'
                },                            
                
            }                           

            ]
        }]
        
        p._remove_redundant_mixing(transfer_instruction_groups)
        
        self.assertEqual(transfer_instruction_groups,
                         [{
                             'transfer': [{
                                 'volume': '900.0:microliter',
                                 'to': 'my_plate/0',
                                 'from': 'culture_medium/0',
                                 'mix_before': {
                                     'volume': '900.0:microliter',
                                     'repetitions': 5,
                                     'speed': '900.0:microliter/second'
                                 }                               
                             },
                             { 
                                 'volume': '900.0:microliter',
                                 'to': 'my_plate/1',
                                 'from': 'culture_medium/0'
                             }  
                                 ]
                         }, {
                             'transfer': [{
                                 'volume': '100.0:microliter',
                                 'to': 'my_plate/0',
                                 'from': 'culture_medium/0'
                             },
                             {
                                 'volume': '900.0:microliter',
                                   'to': 'my_plate/1',
                                   'from': 'culture_medium/1',
                                   'mix_after': {
                                       'volume': '900.0:microliter',
                                       'repetitions': 5,
                                       'speed': '900.0:microliter/second'
                                   },                            
                                   
                               }                                          
                            ]
                         }]
                         )    
        
    def test_remove_redundant_mixing_after(self):
        p = Protocol()
        
        transfer_instruction_groups = [{
                             'transfer': [{
                                 'volume': '900.0:microliter',
                                 'to': 'my_plate/0',
                                 'from': 'culture_medium/0',
                                 'mix_after': {
                                     'volume': '900.0:microliter',
                                     'repetitions': 5,
                                     'speed': '900.0:microliter/second'
                                 }
                             }, {
                                 'volume': '100.0:microliter',
                                 'to': 'my_plate/0',
                                 'from': 'culture_medium/0',
                                 'mix_after': {
                                     'volume': '900.0:microliter',
                                     'repetitions': 5,
                                     'speed': '900.0:microliter/second'
                                 }
                             }]
        }]
        
        p._remove_redundant_mixing(transfer_instruction_groups)
        
        self.assertEqual(transfer_instruction_groups,
                         [{
                             'transfer': [{
                                'volume': '900.0:microliter',
                                'to': 'my_plate/0',
                                'from': 'culture_medium/0'
                               
                                }, {
                                    'volume': '100.0:microliter',
                                    'to': 'my_plate/0',
                                    'from': 'culture_medium/0',
                                    'mix_after': {
                                        'volume': '900.0:microliter',
                                        'repetitions': 5,
                                        'speed': '900.0:microliter/second'
                                    }                                 
                                }]
                         }]
                         )   
        
    def test_remove_redundant_mixing_after_small_volume(self):
            p = Protocol()
            
            transfer_instruction_groups = [{
                                 'transfer': [{
                                     'volume': '9.0:microliter',
                                     'to': 'my_plate/0',
                                     'from': 'culture_medium/0',
                                     'mix_after': {
                                         'volume': '900.0:microliter',
                                         'repetitions': 5,
                                         'speed': '900.0:microliter/second'
                                     }
                                 }, {
                                     'volume': '9.0:microliter',
                                     'to': 'my_plate/0',
                                     'from': 'culture_medium/0',
                                     'mix_after': {
                                         'volume': '900.0:microliter',
                                         'repetitions': 5,
                                         'speed': '900.0:microliter/second'
                                     }
                                 }]
            }]
            
            p._remove_redundant_mixing(transfer_instruction_groups)
            
            #unchanged from original
            self.assertEqual(transfer_instruction_groups,
                             [{
                                 'transfer': [{
                                     'volume': '9.0:microliter',
                                     'to': 'my_plate/0',
                                     'from': 'culture_medium/0',
                                     'mix_after': {
                                         'volume': '900.0:microliter',
                                         'repetitions': 5,
                                         'speed': '900.0:microliter/second'
                                     }
                                 }, {
                                     'volume': '9.0:microliter',
                                     'to': 'my_plate/0',
                                     'from': 'culture_medium/0',
                                     'mix_after': {
                                         'volume': '900.0:microliter',
                                         'repetitions': 5,
                                         'speed': '900.0:microliter/second'
                                     }
                                 }]
            }]
                             )       
        
        
    def test_trash_max_volume(self):
        
        p = Protocol()
        
        plate1 = p.ref("my_plate", cont_type="24-deep", 
                                     storage="cold_4", discard=False)
        
        plate2 = p.ref("my_plate2", cont_type="24-deep", 
                                             storage="cold_4", discard=False)        
        
        for well in plate1.all_wells():
            well.volume = get_well_max_volume(well)
           
        p.transfer(plate1.well(0),plate2.well(0),ul(100))
        
        p.trash_max_volume(plate1.all_wells())
        
        assert all([well.volume==get_well_dead_volume(well) for well in plate1.all_wells()])
        
        self.assertEqual(len(p.instructions),2)
        
    def test_consolidate_transfer(self):
        
        p = Protocol()
        
        plate1 = p.ref("my_plate", cont_type="24-deep", 
                                     storage="cold_4", discard=False)
        
        plate2 = p.ref("my_plate2", cont_type="24-deep", 
                                             storage="cold_4", discard=False)        
        
        for well in plate1.all_wells():
            well.volume = get_well_max_volume(well)
        
        
           
        p.transfer(plate1.well(0),plate2.well(0),ul(100))
        p.transfer(plate1.well(0),plate2.well(0),ul(100))
        p.distribute(plate1.well(1), [plate2.well(1)], ul(100))
        p.transfer(plate1.well(0),plate2.well(0),ul(100))
        
        p._consolidate_last_pipette_operation()
        
        self.assertListEqual(p.as_dict(seal_on_store=False)['instructions'],
                             [{
                                 'groups': [{
                                     'transfer': [{
                                         'volume': '100.0:microliter',
                                         'to': 'my_plate2/0',
                                         'from': 'my_plate/0'
                                         }, {
                                             'volume': '100.0:microliter',
                                             'to': 'my_plate2/0',
                                             'from': 'my_plate/0'
                                         }]
                                     }, {
                                         'distribute': {
                                             'to': [{
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/1'
                                                 }],
                                             'from': 'my_plate/1',
                                             'allow_carryover': False
                                         }
                                         }, {
                                             'transfer': [{
                                                 'volume': '100.0:microliter',
                                                 'to': 'my_plate2/0',
                                                 'from': 'my_plate/0'
                                             }]
                                             }],
                                 'op': 'pipette'
                             }]
                             )
        
    def test_transfer_allows_containers(self):
            
        p = Protocol()
        
        tube = p.ref("my_tube", cont_type="micro-2.0", 
                                     storage="cold_4", discard=False)
        
        plate2 = p.ref("my_plate2", cont_type="6-flat", 
                       storage="cold_4", discard=False)        
        
        for well in tube.all_wells():
            well.volume = get_well_max_volume(well)
           
        p.transfer(tube,plate2,ul(100))
        
        
        assert all([well.volume==ul(100) for well in plate2.all_wells()])
        
        self.assertEqual(tube.well(0).volume,ul(1400))    
    
        
    def test_distribute(self):
            
            p = Protocol()
            
            plate1 = p.ref("my_plate", cont_type="24-deep", 
                                         storage="cold_4", discard=False)
            
            plate2 = p.ref("my_plate2", cont_type="24-deep", 
                                                 storage="cold_4", discard=False)        
            
            for well in plate1.all_wells():
                well.volume = get_well_max_volume(well)
               
            p.distribute(plate1.well(0),plate2.all_wells(),ul(100),
                         mix_before=True
                         )
            
            self.assertEqual(p.as_dict(seal_on_store=False),
                             {
                                 'refs': {
                                     'my_plate2': {
                                         'new': '24-deep',
                                         'store': {
                                             'where': 'cold_4'
                                         }
                                     },
                                     'my_plate': {
                                         'new': '24-deep',
                                         'store': {
                                             'where': 'cold_4'
                                         }
                                     }
                                 },
                                 'instructions': [{
                                     'groups': [{
                                         'distribute': {
                                             'to': [{
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/0'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/1'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/2'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/3'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/4'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/5'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/6'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/7'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/8'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/9'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/10'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/11'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/12'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/13'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/14'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/15'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/16'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/17'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/18'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/19'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/20'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/21'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/22'
                                             }, {
                                                 'volume': '100.0:microliter',
                                                 'well': 'my_plate2/23'
                                             }],
                                             'from': 'my_plate/0',
                                             'mix_before': {
                                                 'volume': '900.0:microliter',
                                                 'repetitions': 5,
                                                 'speed': '900.0:microliter/second'
                                             },
                                             'allow_carryover': False
                                         }
                                     }],
                                     'op': 'pipette'
                                 }]
                             }
                             )
        
        
    
    def test_provision_by_name_allows_containers(self):
        
        initial_protocol = Protocol()
        
        plate = initial_protocol.ref("my_plate", cont_type="24-deep", 
                                     storage="cold_4", discard=False)
        
        initial_protocol.provision_by_name('water', plate, ml(1))
        
        self.assertTrue( all([well.volume==ml(1) for well in plate.all_wells()]) )
        

    def test_provision_by_name_transcriptic_reagent_warming(self):
        
        initial_protocol = Protocol(mammalian_cell_mode=True)
        
        plate = initial_protocol.ref("my_plate", cont_type="24-deep", 
                                     storage="cold_4", discard=False)
        
        initial_protocol.provision_by_name(Reagent.trypsin_edta, plate.well(0), ml(1),
                                           pre_warm_minutes=5)
        
        p = initial_protocol._get_final_protocol()
        
        #check that the instructions include provisioning statements at the beginning for
        assert isinstance(p.instructions[0],Provision)
        
        self.assertEqual([instruction.data['resource_id'] for instruction in p.instructions[:1]],
                         [p.transcriptic_inv[Reagent.trypsin_edta]])

            
        #final protocol needs to be the same as the initial protocol
        self.assertEqual(p.as_dict(seal_on_store=False),initial_protocol.as_dict(seal_on_store=False))
        self.assertDictEqual(p.as_dict(seal_on_store=False),
                             {
                                 'refs': {
                                     'trypsin_edta': {
                                         'new': '96-deep',
                                         'discard': True
                                     },
                                     'my_plate': {
                                         'new': '24-deep',
                                         'store': {
                                             'where': 'cold_4'
                                         }
                                     }
                                 },
                                 'time_constraints': [{
                                     'to': {
                                         'instruction_start': 3
                                     },
                                     'less_than': '1.0:hour',
                                     'from': {
                                         'instruction_end': 1
                                     }
                                 }, {
                                     'to': {
                                         'instruction_end': 1
                                     },
                                     'less_than': '1.0:hour',
                                     'from': {
                                         'instruction_start': 3
                                     }
                                 }],
                                 'instructions': [{
                                     'to': [{
                                         'volume': '900.0:microliter',
                                         'well': 'trypsin_edta/0'
                                     }, {
                                         'volume': '115.0:microliter',
                                         'well': 'trypsin_edta/0'
                                     }],
                                     'op': 'provision',
                                     'resource_id': 'rs196bb2deuhsy'
                                 }, {
                                     'lid': 'standard',
                                     'object': 'trypsin_edta',
                                     'op': 'cover'
                                 }, {
                                     'where': 'warm_37',
                                     'object': 'trypsin_edta',
                                     'co2_percent': 0,
                                     'duration': '5:minute',
                                     'shaking': False,
                                     'op': 'incubate'
                                 }, {
                                     'object': 'trypsin_edta',
                                     'op': 'uncover'
                                 }, {
                                     'groups': [{
                                         'transfer': [{
                                             'volume': '900.0:microliter',
                                             'to': 'my_plate/0',
                                             'from': 'trypsin_edta/0'
                                         }, {
                                             'volume': '100.0:microliter',
                                             'to': 'my_plate/0',
                                             'from': 'trypsin_edta/0'
                                         }]
                                     }],
                                     'op': 'pipette'
                                 }]
                             }                        
                         )
        
        
    def test_provision_by_name_adds_time_constraints(self):
            
        p = Protocol()
        
        plate = p.ref("my_plate", cont_type="24-deep", 
                                     storage="cold_4", discard=False)
        
        p.provision_by_name('culture_medium', plate.well(0), ml(1),
                            pre_warm_minutes=5)
        
        incubate_instructions = [instruction for instruction in p.instructions \
                                 if isinstance(instruction, Incubate)]   
        
        self.assertEqual(len(incubate_instructions),1)
        
        self.assertEqual(incubate_instructions[0].data['duration'],
                         '5:minute')
        
        
        incubate_instruction_index = None
        
        for i in range(0,len(p.instructions)):
            if p.instructions[i] == incubate_instructions[0]:
                incubate_instruction_index = i
                break
            
        pipette_instruction_index = None
                
        for i in range(incubate_instruction_index,len(p.instructions)):
            if p.instructions[i].op == 'pipette':
                pipette_instruction_index = i
                break        
            
        self.assertEqual(len(p.time_constraints),2)
        self.assertEqual(p.time_constraints,[{
                  "to": {
                    "instruction_start": pipette_instruction_index
                  },
                  "less_than": hours(1),
                  "from": {
                    "instruction_end": incubate_instruction_index
                  }
                },
        
                {
                    "to": {
                        "instruction_end": incubate_instruction_index
                    },
                    "less_than": hours(1),
                    "from": {
                        "instruction_start": pipette_instruction_index
                    }
                  }
                ])  
        
    def test_transfer_all_volume_evenly(self):
        
        p = Protocol()
        
        source_plate = p.ref("source_plate", cont_type="6-flat", 
                                            storage="cold_4", discard=False)
        
        for well in source_plate.all_wells():
            well.volume = ml(3)
        
        dest_tubes = []
        for i in range(0,8):
            dest_tubes.append(p.ref('freeze_tube_%s'%i, cont_type='micro-2.0', 
                                     discard=True))        
            
        p.transfer_all_volume_evenly(source_plate, 
                                    dest_tubes)
        
        self.assertTrue( all([tube.well(0).volume==ul(1950) for tube in dest_tubes]) )
        self.assertTrue( all([well.volume==ul(400) for well in source_plate.all_wells()]) )
        
        
    def test_ref_custom_property(self):

        p = Protocol()
        
        plate = p.ref("plate", cont_type="6-flat", 
                      storage="cold_4", discard=False,
                      cell_line_name = 'vero_1',
                      custom_property2 = 'test'
                      )        
        self.assertTrue( all([well.properties['cell_line_name']=='vero_1' for well in plate.all_wells()]) )
        self.assertTrue( all([well.properties['custom_property2']=='test' for well in plate.all_wells()]) )
            
        
    def test_spin(self):
        
        p = Protocol()
        
        cell_plate = create_blank_plate('96-flat')
        
        for well in cell_plate.all_wells():
            well.volume = get_well_max_volume(well)        
        
        p.spin(cell_plate, '8.4:g', '2:second',flow_direction='outward')
        
        self.assertTrue(all([well.volume == get_well_dead_volume(well) for well in cell_plate.all_wells()]))
    
    def test_add_antibiotic(self):
    
        p = Protocol()
    
        plate = create_blank_plate('96-deep')
        well = plate.well(0)
    
        well.volume = ml(1)
    
        p.add_antibiotic(well, Antibiotic.amp)
    
        self.assertEqual(well.volume, ul(1001))
    
    
        well = plate.well(1)
    
        well.volume = ul(615)
    
        p.add_antibiotic(well, Antibiotic.amp)
    
        self.assertEqual(well.volume, ul(615.7))    
        
        
        well = plate.well(2)
    
        well.volume = ul(615)
    
        p.add_antibiotic(well, Antibiotic.strep)
    
        self.assertEqual(well.volume, ul(615.7))         
        
    def test_measure_bacterial_density(self):
        p = Protocol()
    
        plate = create_blank_plate('96-deep')
        wells = get_column_wells(plate,0)
        for well in wells:
            well.volume=ul(200)
            
        p.measure_bacterial_density(wells)
        
        self.assertEqual(type(p.instructions[p.get_instruction_index()]),Absorbance)



        
    def test_complete_mix_kwargs(self):
        
        p = Protocol()
                
        cell_plate = create_blank_plate('96-flat')    
        
        
        cell_plate.well(0).volume = ul(340)
        cell_plate.well(1).volume = ul(0)
          
          
        #test basic  
        mix_kwargs = {
            'mix_after': True,
            'mix_before': True
        }          
        
        p._complete_mix_kwargs(mix_kwargs, cell_plate.wells(0), cell_plate.wells(1), ul(100),
                               allow_different_mix_options=False)
        
        self.assertDictEqual(mix_kwargs,
                             {
                                 'mix_vol': ul(50),
                                 'mix_before': True,
                                 'flowrate': '50.0:microliter/second',
                                 'mix_after': True,
                                 'repetitions': 5
                             }                             
                             )        
        
        #test diff mix options
        mix_kwargs = {
                    'mix_after': True,
                    'mix_before': True
                }          
        
        cell_plate.well(0).volume = ul(340)
        cell_plate.well(1).volume = ul(0)        
                
        p._complete_mix_kwargs(mix_kwargs, cell_plate.wells(0), cell_plate.wells(1), ul(100),
                               allow_different_mix_options=True)
        
        self.assertDictEqual(mix_kwargs,
                             {
                                 'mix_vol_a': ul(50),
                                 'mix_vol_b': ul(170),
                                 'flowrate_a': '50.0:microliter/second',
                                 'flowrate_b': '170.0:microliter/second',
                                 'repetitions_a': 5,
                                 'repetitions_b': 5,
                                 'mix_before': True,
                                 'mix_after': True
                             }                             
                             )   
        
        #test mix_seconds
        mix_kwargs = {
                    'mix_before': True,
                    'mix_seconds': 4,
                    'mix_percent': 25
                }          
        
        cell_plate.well(0).volume = ul(340)
        cell_plate.well(1).volume = ul(0)  
        
        p._complete_mix_kwargs(mix_kwargs, cell_plate.wells(0), cell_plate.wells(1), ul(100),
                               allow_different_mix_options=False)
        
        p.transfer(cell_plate.wells(0), cell_plate.wells(1), ul(100), mix_seconds=4, mix_after=True)
        
        self.assertDictEqual(mix_kwargs,
                             {
                                 'mix_vol': ul(85),
                                 'flowrate': '42.0:microliter/second',
                                 'mix_before': True,
                                 'repetitions': 5
                             }                             
                             )  
        #test mix_seconds
        mix_kwargs = {
                    'mix_after': True,
                    'mix_before': True,
                    'mix_seconds_b': 4,
                    'mix_percent_b': 75,
                    'repetitions_b': 2,
                    'mix_seconds_a': 4,
                    'repetitions_a': 3,
                    'mix_percent_a': 50                    
                }          
        
        cell_plate.well(0).volume = ul(340)
        cell_plate.well(1).volume = ul(0)        
                
        p._complete_mix_kwargs(mix_kwargs, cell_plate.wells(0), cell_plate.wells(1), ul(100),
                               allow_different_mix_options=True)
        
        self.assertDictEqual(mix_kwargs,
                             {
                                 'mix_vol_a': ul(50),
                                 'mix_vol_b': ul(255),
                                 'flowrate_a': '25.0:microliter/second',
                                 'flowrate_b': '127.0:microliter/second',
                                 'repetitions_a': 3,
                                 'repetitions_b': 2,
                                 'mix_before': True,
                                 'mix_after': True
                             }                             
                             )         
       
        
        


        


    def test_dispense_by_name(self):
        p = Protocol()
        plate1 = p.ref('test plate', cont_type='24-deep', discard=True)
        
        p.dispense_by_name(Reagent.lb_amp, plate1, [{"column": 0, "volume": ml(2)},
            {"column": 0, "volume": ml(1)},
            {"column": 1, "volume": ml(2)},
            {"column": 2, "volume": ml(1)}])
        
        self.assertTrue(all([well.volume == ml(3) for well in get_column_wells(plate1,0)]))
        self.assertTrue(all([well.volume == ml(2) for well in get_column_wells(plate1,1)]))
        self.assertTrue(all([well.volume == ml(1) for well in get_column_wells(plate1,2)]))
        self.assertTrue(all([well.volume == ml(0) for well in get_column_wells(plate1,3)]))
    
    def test_dispense_by_name_volume_check(self):
        p = Protocol()
        plate1 = p.ref('test plate', cont_type='96-pcr', discard=True)
        
        
        with self.assertRaises(Exception):
            p.dispense_by_name(Reagent.lb_amp, plate1, [{"column": 0, "volume": ml(2)}])
    
    def test_stamp_volume_check(self):
        pass
    
    def test_transfer_column_single_column(self):
        p = Protocol()
        plate1 = p.ref('test plate', cont_type='96-pcr', discard=True)        
        plate2 = p.ref('test plate 2', cont_type='96-pcr', discard=True)
        
        for well in plate1.all_wells():
            well.volume = ul(160)
            
        p.transfer_column(plate1, 0, 
                         plate2, 
                         1, 
                         ul(50), mix_before=True)
        
        self.assertTrue(all([well.volume==ul(50) for well in get_column_wells(plate2,1)]))
        self.assertTrue(all([well.volume==ul(110) for well in get_column_wells(plate1,0)]))
        
    
    def test_transfer_column_multi_column(self):
        p = Protocol()
        plate1 = p.ref('test plate', cont_type='96-pcr', discard=True)        
        plate2 = p.ref('test plate 2', cont_type='96-pcr', discard=True)
    
        for well in plate1.all_wells():
            well.volume = ul(160)
    
        p.transfer_column(plate1, 0, 
                          plate2, 
                          [0,1,1,2], 
                          ul(20), mix_before=True,
                          one_tip=True)
    
        self.assertTrue(all([well.volume==ul(20) for well in get_column_wells(plate2,0)]))
        self.assertTrue(all([well.volume==ul(40) for well in get_column_wells(plate2,1)]))
        self.assertTrue(all([well.volume==ul(20) for well in get_column_wells(plate2,0)]))
        self.assertTrue(all([well.volume==ul(80) for well in get_column_wells(plate1,0)]))    
        
        p.transfer_column(plate1, [1,2], 
                          plate2, 
                          [3,4], 
                          ul(65), mix_before=True,
                          one_tip=True)        
        
        self.assertTrue(all([well.volume==ul(65) for well in get_column_wells(plate2,[3,4])]))
        self.assertTrue(all([well.volume==ul(95) for well in get_column_wells(plate1,[1,2])])) 
        
        
    def test_stamp_zero_volume(self):
        p = Protocol()
        