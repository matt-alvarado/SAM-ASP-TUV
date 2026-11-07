#!/usr/bin/env python3

MAX_STRING_LENGTH = 16

def read_file(f):
    """ TODO """
    list_of_species = []
    try:
        with open(f, 'r') as fptr:
            for line in fptr:
                if (line[0] != '!' and line[0] != '0' and line[0] != '1'
                        and "END-OF-FILE" not in line): 
                    temp_line = line.strip().split(';')
                    list_of_species.append(temp_line[0])
    except IOError:
        raise
    finally:
        return list_of_species

def write_template(list_of_species, outfilename, newline='\n'):
    """ todo """
    try:
        with open(outfilename, 'w') as fptr:
            for species in list_of_species:
                fptr.write("%16s;%s" % (species.strip(), newline))
    except IOError:
        raise

species = read_file("../GasPhaseChems.in")
write_template(species, "_GaseousSpeciesBlankList.in")
        
    
