#!/usr/bin/env python3

import numpy as np
import sys
import csv 

num_rxns = 51
sub_tuple = (num_rxns,)
#sub_str = '%df8' % num_rxns
zenith_angle_bins = [[0, 10], [10, 20], [20, 30], [30, 40], [40, 50],
                     [50, 60], [60, 70], [70, 76], [76, 86], [86, 90]]
dt = [('utc_hour', float), ('zenith_angle', float), 
      ('reactions', np.float64, sub_tuple)]#('%d,'% num_rxns))]
print(dt)
def read_file_data(f='/rotor/scratch/cbrodows/ASP_STILT-Chem_Mapping/InorganicInputDecks/PhotoRates.in'):
    """ TODO """
    list_of_tuples = []
    try:
        with open(f) as fptr:
            for line in fptr:
                if line.startswith('!') or line.startswith('END'):
                    continue
                else:
                    wanted_line = line.strip()[:-1]
                    parts = wanted_line.split()
                    new_parts = []
                    new_parts.append(parts.pop(0)); new_parts.append(parts.pop(0))
                    new_parts.append(tuple(parts))
                    #print 'new_parts = ', new_parts
                    #print 'length = ', len(new_parts)
                    #print 'length last section = ', len(new_parts[2])
                    #sys.exit()
                    list_of_tuples.append(tuple(new_parts))
        #print 'list of tuples = ', list_of_tuples
        #sys.exit()
        return np.array(list_of_tuples, dtype=dt)
        #    lines = [line for line in fptr if not line.startswith('!') and not line.startswith('END'))
        #    arr = np.genfromtxt(lines, comments='!', dtype=dt)
        #return arr
    except:
        raise
        sys.exit()

def get_species(f, header_lines=16):
    """
    """
    linectr = 0
    species_dict = dict()
    try:
        with open(f, 'r') as fptr:
            
            ln = fptr.readline()
            while 'Photolysis' not in ln:
            #while linectr < header_lines:
                #print 'Skipping header...'
                ln = fptr.readline()
                linectr += 1
            print('past header...last line = ', ln)
            line = fptr.readline()  
            #if 'Photolysis' in line2:
            #    #print 'Skipping Photolysis rate line'
            #    line = fptr.readline()
            #    linectr += 1
            
            #line = fptr.readline() 
            while True:
            #print 'fptr.readline = ', fptr.readline() 
            #for line in fptr.readlines():
                #print '\nLine = ', line
                line_without_ex = line[1:].strip()
                #print 'Line = ', line_without_ex
                if "Values" in line_without_ex or 'values' in line_without_ex:
                    values = False
                    break

                fields = line_without_ex.split()
                print(fields)                
                equation_num = fields[0]
                species_name = fields[2]
                
                if species_name in species_dict:
                    species_copy = species_name + "<2>"
                    clone_num = 2
                     
                    while species_copy in species_dict:
                        clone_num += 1
                        species_copy = species_name + "<%d>" % clone_num
                    
                    # Assign back for further processing    
                    species_name = species_copy
                
                # Store the values as strings so that the dictionary can be
                # easily inverted    
                species_dict[species_name] = str(equation_num)
                
                line = fptr.readline()
                if 'END-' in line:
                    break
            
    except IOError:
        raise
    finally:
        return species_dict
        
def read_file(f=('/rotor/scratch/cbrodows/ASP_STILT-Chem_Mapping' +
                 '/InorganicInputDecks/PhotoRates.in')):
    """ TODO
    """
    array = read_file_data(f)
    species_dict = get_species(f)
    #print species_dict; sys.exit()
    return array, species_dict
    
def bin_data_by_zenith_angle(data, species_dict, result='mean'):
    """
    Current Implementation takes all the reactions at a given zenith angle 
    bin and will take the mean (or something.)

    TODO: Make as rations w/NO2...
    """
    lines = []
    #no2_factor = 1.0; no2_flag = False
    lines.append(' %d\n' % num_rxns)
    #inverse_dict = dict(zip(map.values(),map.keys()))
    #inverse_dict = {value: key for key, value in species_dict.items()}
    inverse_dict = dict((v,k) for k, v in species_dict.items()) 
    #{value: key for key, value in species_dict.items()}
    print([key for key in inverse_dict])
    no2_index = species_dict['NO2']
    reactions = ([no2_index]+list([species_dict[key] 
                  for key in species_dict if key != 'NO2']))
    no2_factors = np.array([1.0 for i in range(len(zenith_angle_bins))])
    for rxn in reactions: #xrange(1, num_rxns + 1, 1):
        rxn_num = int(rxn)
        species = inverse_dict[str(rxn_num)]
        print('Reaction #: ', rxn_num)
        print('Species: ', species)

        #if not no2_flag:

        # Get the flag
        flag = 1
        if '<' in species:
            parts1 = species.strip().split('<')
            species_root = parts1[0]
            parts2 = parts1[-1].strip().split('>')
            flag = int(parts2[0])
                    
        string = " %s %d" % (species, flag) #inverse_dict[str(rxn_num)])
        
        zenith_bin_ctr = -1
        for zenith_bin in zenith_angle_bins:
            zenith_bin_ctr += 1
            print('\tZenith Bin: %d to %d' % (zenith_bin[0], zenith_bin[1]))
            bin_data = data[np.where((data['zenith_angle'] >= 
                                      zenith_bin[0]) &
                                     (data['zenith_angle'] < 
                                      zenith_bin[1]))]
            #no2_factor = 1.0
            if len(bin_data) == 0:
                rxn_constant = 0.0#-999.9
            else:
                    try:
                            #print '\tAvailable Data: ', bin_data['reactions']
                            #sys.exit()
                            rxn_constants = np.array([arr[rxn_num - 1] for arr in bin_data['reactions']])
                            print('\tAvailable Data: ', rxn_constants)
                        
                            if isinstance(rxn_constants, np.ndarray):
                                if result == 'mean':
                                    rxn_constant = np.mean(rxn_constants)
                                elif result == 'median':
                                    rxn_constant = np.median(rxn_constants)
                                    
                                # Interpolate?
                                else:
                                    assert False, ("TODO: mean/median only numpy " +
                                                   "routines so far!")
                            else:
                                rxn_constant = rxn_constants

                            if species == 'NO2':
                                no2_factors[zenith_bin_ctr] = rxn_constant
                                #no2_factor = rxn_constant

                    except IndexError:
                            print('no2_factors = ', no2_factors)
                            print('IndexError; rxn num = ', rxn_num)
                            #rxn_constant = 0.0#-999.9
                    
            #string = string + " %.5f" % (rxn_constant/no2_factor)
            print('prior to string write: no2_factors = ', no2_factors)
            print('desired value: no2_factors[%d] = %f' % (zenith_bin_ctr, no2_factors[zenith_bin_ctr]))
            string = string + " %.3E" % (rxn_constant/no2_factors[zenith_bin_ctr])
        string = string + "\n"
        lines.append(string)
    
    # Let's add a commented line containing the NO2 factors
    str_temp = "! NO2 totals:"
    for no2_factor in no2_factors:
        str_temp += ' %.3E' % no2_factor
    str_temp += '\n'
    lines.append(str_temp)
    return lines
    
def write_to_file(lines, outf="PHOTO.DAT"):
    """ TODO
    """
    try:
        with open(outf, 'w') as fptr:
            fptr.writelines(lines)
    except IOError:
        raise
        
array, d = read_file("../PhotoRates.generic.200DU.1km.in")        
write_to_file(bin_data_by_zenith_angle(array, d))
sys.exit()                
            
                
    
