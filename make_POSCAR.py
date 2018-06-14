import numpy as np
from operator import mul
from sklearn.cluster import KMeans
import math

####################################################################
######################## ISSUES START ##############################
####################################################################
# 1 - Need to make period of changes more robust
# 2 - Make function to automatically get all the material details
# 3 - Pack everything into classes
# 4 - Stick to data abstraction
# 5 - Fill in doctests everywhere
# 6 - Make stacking types more general?
# 7 - add more edge types to ribbon edge methods, improve the method
# 8 - cleanup unit cell maker adustment for chopped ribbons

# 9 - fix the stacking offset - it works but its mixed up AC and AD for some reason - something to do with cycle?
      # fudged it for now and works if you define AA etc at the bottom
# 10 - Lithium intercalation
####################################################################
######################### ISSUES END ###############################
####################################################################





####################################################################
################# DATA ABSTRACTION START ###########################
####################################################################

def pick_column(array,axis_index):
    """
    This is a selector to pick the axes of an atom data structure
    :param array: a numpy array
    :param axis_index: an integer which is used to select an axis of the array
    :return: a numpy array of size n in 0th dimension and 1 in 1st dimension
    """

    return array[:,axis_index]


def offset_builder(stacking, atoms_per_cell, yoffset):
    """
    This builds a structure of offsets in each direction given the stacking type
    :param stacking: a string to indicate the type of stacking, e.g. 'AB'
    :return: a 1x3 numpy array with each column as one of the axis offsets in terms of unit cells.
    This is the shift in a
    monolayer upon stacking, e.g. a bilayer is one monolayer plus one monolayer offset by these coordinates.
    We do it 4 times to account for each of the four p atoms in unit cell
    """

    if stacking == 'AA':
        # offset of layers #This is defined from converged bulk pe-d2 vdw
        xoffset = 0
        zoffset = 0
    elif stacking == 'AB':
        xoffset = 0.5
        zoffset = 0
    elif stacking == 'AC':
        xoffset = 0.5
        zoffset = 0.5
    elif stacking == 'AD':
        xoffset = 0
        zoffset = 0.5
    else:
        raise ValueError("please enter a valid stacking type. Only supported stackings are AA,AB,AC,AD")

    return repeater( array_builder(xoffset, yoffset, zoffset), atoms_per_cell)


def array_builder(i,j,k):
    """
    A function to return a numpy array 1x3 dimension with i, j, k as elements.
    :return: np array
    """

    return np.array([[i, j, k]])

def stacker(unit, ln):
    """
    A function to build ln stacks of the same array
    :param unit: the unit you wish to stack
    :param ln: the number of layers you wish to stack
    :return: an array which is ln repeats of unit
    """

    return np.tile(unit,(ln,1))

def repeater(unit,thickness):
    """
    A function which repeats layers - not stacking. e.g 1 2 -> 1 1 2 2 (if these are vertical arrays)
    :param unit: the unit which will be repeated
    :param thickness: the thickness of the repeat. e.g. 2 in the example above
    :return: a np array of thickness repeated units
    """

    return np.repeat(unit,thickness,axis=0)



def cycle(period,length):
    """
    ###ABSTRACTION BARRIER!!!###

    A function which builds a np vector which cycles between 0 and 1 depending on the period of repition
    :param period: a number corresponding periodicity of 0s
    :param length: the lnegth of the vector which should be output
    :return: an np vector of length l repeating 0s every period and 1s elsewhere
    """


    cycled = np.zeros((length,1))
    for i in range(length):
        cycled[i] = bool((i%period))

    return cycled # A list of 1's and 0's of period period.


def shifter(array,period,atoms_per_cell,ln,axes):
    """
    ###This breaks data abstraction barrier because you can't assign a function to a function!!###\

    This is a function which builds a shifted array. This takes in an array which is stacked, and shifts each layer
    by an appropriate amount depending on the period given. i.e. cycles periodic ones or scales stacking offsets
    :param array: a stacked np array of the offsets
    :param period: the period np array of period of shifts
    :param atoms_per_cell: the number of atoms in a unit cell, i.e. how many lines of the array belong to each monolayer
    :param ln: the number of layers we will stack
    :return: a np array which is appropriately shifted, but not wrapped
    """

    for i in range(axes):

        if pick_column(period,i): # if that column is a number then adjust accordingly
            multiplier = repeater(cycle(period[0, i], ln), atoms_per_cell) # alternates based on period
            array[:,i] = np.multiply(array[:,i],multiplier[:,0]) # replaces ith column with the offset multiplied by multiplier

        else:
            multiplier = repeater(np.arange(ln).reshape(ln, 1),atoms_per_cell)[:,0] # scales
            array[:, i] = np.multiply(array[:, i], multiplier) # replaces ith column with the offset multiplied by multiplier

    return array #returns an array with columns shifted by necessary amounts/cycled if they have a period


def wrapper(atoms,period,lattice_const):
    """
    ####ABSTRACTION BARRIER BROKEN####
    A function which wraps around shifts which are supposed to be periodic. e.g. shifting x by 0.5 each time
    should lead to wrapping 1 unit cell shift back to zero

    :param atoms: a np array of shifted values which now needs to be wrapped
    :param period: the period array of shifts
    :param lattice_const: the lattice constant array
    :return: a np array which is now properly wrapped
    """

    for i in range(axes):
        if period[0,i]:
            atoms[:, i] = (atoms[:,i]%1)*lattice_const[0,i]

    return atoms


def adder(data_struct1, data_struct2):
    """
    A data abstraction ot combine two atom structures. This currently assumes they are np arrays
    :param data_struct1: an np array
    :param data_struct2: an np array
    :return: an np array
    """

    return data_struct1 + data_struct2


def combine(struct1, struct2):
    """
    multiply together two np arrays element-wise
    :param struct1: np array
    :param struct2: np array
    :return: np array of two inputs multiplied together element wise
    """

    combined = struct1 * struct2

    return combined


def permuter(xw=1, zw=1):
    """
    This loop creates an array with permutations of xw and zw unit cell positions throughout.
    :param xw: the number of permutations in the first column
    :param zw: the number of permutations in the second column
    :return: an np array with all combinations e.g. [0, 0, 0] [1, 0, 0] etc
    """
    base = np.array([[]])
    for i in range(xw):
        j = 0
        basex = np.array([[i, 0, j]])
        try:
            base = np.append(base, basex, axis=0)
        except:
            base = np.array([[0, 0, 0]])
        for j in range(zw - 1):
            basez = np.array([[i, 0, j + 1]])
            base = np.append(base, basez, axis=0)

    return base


def add_buffer(array_col, spacer, width):
    """
    Adds a buffer to atom positions or unit cell only if the width is greater than one unit cell
    :param array_col: the column of the array to be added to
    :param spacer: the size of the spacer to be added
    :param width: the
    :return:
    """

    buffed = array_col + spacer * bool(width-1) # add only if more than one unit cell

    return buffed


####################################################################
################# DATA ABSTRACTION END #############################
####################################################################







####################################################################
####################### FUNC DETAILS START #########################
####################################################################


def monolayer_builder(spacer, thickness):
    """
    A funciton to build an np array of monolayer atom positions
    :return: np array of monolayer atom positions
    """

    #defined from converged bulk pe-d2 vdw

    x1 = 0
    x2 = 0.500000
    y1 = spacer
    y2 = spacer + thickness
    z1 = 0.918924
    z2 = 0.581076
    z3 = 0.081076
    z4 = 0.418924

    mono = np.array([[x1, y1, z1], [x2, y1, z2], [x1, y2, z3], [x2, y2, z4]])

    return mono



def atom_stacker(mono, ln, atoms_per_cell, axes, lattice_const, yoffset, stacking, pera, perb, perc):
    """
	A function which takes in a monolayer and number of layers, and outputs a np array which is shifted correctly for each layer
	:param mono: np array of monolayer base unit
	:param ln: integer of number of layers
	:return: np array of all atom positions
	"""

    offset = offset_builder(stacking, atoms_per_cell, yoffset) # This returns a data structure with the offset information for each layer in.
    period = array_builder(pera, perb, perc) # Return a data structure with information on the period of an offset.
                              # if no period then 0 (e.g. y axis stacking). If e.g. x axis offset period is every other
                              # monolayer then period would be 2
                             # 2 0 0  means x repeats every 2 cells e.g.ABABA

    repeated_offset = stacker(offset, ln)  # This returns a data structure which is ln repeats of offset
    shifted_offset = shifter(repeated_offset, period, atoms_per_cell,ln, axes) # Return the offset with the period of repeats accounted for
    repeated_monolayers = stacker(mono, ln)  # This returns ln repeats of monolayer.
    final_atoms = adder(repeated_monolayers, shifted_offset) # Shift the repeated monolayers by the appropriate amount
    final_atoms = wrapper(final_atoms, period, lattice_const)

    return final_atoms




#use abstraction functions
def cell_expander(a,b,c,ln,ygap,spacing,xw=1,zw=1,distance_chopped=0):
    """
    Function which expands the cell size correctly.
    :param a: a parameter in A
    :param b: layer thickness
    :param c: c parameter in A
    :param ln: integer number of layers
    :param xw: width of the ribbon in unit cells. 1=no ribbon created
    :param zw: width of the ribbon in unit cells. 1=no ribbon created
    :return: a 3x3 np array of cell dimensions for use in POSCAR format for conventional cell
    """

    x_width = a*xw + 2*spacing*bool(xw - 1) # unit cells width plus padding only if making a ribbon
    y_height = 10 + 10 + ygap*(ln - 1) + (thickness)*ln  # bottom padding + top padding+ total gaps between layers + total thickness of layers
    z_width = c*zw + 2*spacing*bool(zw - 1) # unit cells width plus padding only if making a ribbon


    #this gets rid of any extra padding leftover from chopping edge atoms
    if (xw-1):
        x_width-=distance_chopped
    elif (zw-1):
        z_width-=distance_chopped



    return np.array([[x_width,0,0],[0,y_height,0],[0,0,z_width]])



#function to build the poscar string
def build_poscar(atoms, unit_cell, ln, atom_type, atoms_per_cell, combo_rib_dim, no_deleted , edge_atom='', edge_atom_no='', edge_atom_positions=np.array([]) ,intercalant_atom='', intercalant_atom_no='', intercalant_atom_positions=np.array([])):
    """
    Fucntion to write the details of atom type, cell dimensions and atom positions to POSCAR format string.
    UPDATE TO be able to write other types of atoms e.g. for Li intercalation? - or make new function which can call this
    :param atoms: atoms position np array in A
    :param unit_cell: cell dimensions np array in A
    :param ln: layer numbers
    :param atom_type: type of atom
    :param atoms_per_cell: number of atoms per unit cell
    :param combo_rib_dim: this is an integer for combination of ribbon widths to calculate the number of atoms in the cell
    :return: a string in POSCAR format
    """

    start_string=(str(atom_type)+' '+str(intercalant_atom)+ ' '+str(edge_atom)+'\n'+"1.0"+'\n') # the initial info in poscar
    mid_string=('\n'+'\t'+str(atom_type)+ ' '+ str(intercalant_atom) + ' ' + str(edge_atom)+'\n'+'\t'+
        str(atoms_per_cell*ln*combo_rib_dim-no_deleted)+'\t'*bool(intercalant_atom_no)+str(intercalant_atom_no)+ '\t'*bool(edge_atom_no)+str(edge_atom_no)+'\n'+"Cartesian"+'\n') # the middle part with some  information
    spaces=" " * 9 # spacing to make it look pretty

    # turn the numpy arrays for unit cell and positions into strings
    unit_cell_string='\n\t\t'.join(spaces.join('%0.10f' %x for x in y) for y in unit_cell)
    cartesian_string='\n\t'.join(spaces.join('%0.10f' %x for x in y) for y in atoms)


    #if intercalant atom array is not empty then turn into a string, otherwise set to blank
    if intercalant_atom_positions.size:
        intercalant_atom_string = '\n\t'.join(spaces.join('%0.10f' %x for x in y) for y in intercalant_atom_positions)
    else:
        intercalant_atom_string = ''


    if edge_atom_positions.size:
        edge_atoms_string = '\n\t'.join(spaces.join('%0.10f' %x for x in y) for y in edge_atom_positions)
    else:
        edge_atoms_string = ''

    print(edge_atoms_string)

    return start_string + '\t\t' + unit_cell_string + mid_string + '\t' + cartesian_string + '\n\t'*bool(intercalant_atom_string) + intercalant_atom_string + '\n\t'*bool(edge_atoms_string)+edge_atoms_string# combine it all



def write_poscar(poscar_string, filename):
    """
    Function to write the POSCAR format string to a file named filename
    :param poscar_string: A POSCAR format string to be written
    :param filename: the filename to be written to
    :return: True
    """

    assert type(filename)==str, "please enter a string as filename"

    filename=filename+".vasp"
    with open(filename,'w') as vasp_file:
        vasp_file.write(poscar_string)

    print("writing completed")
    return True


def ribbon_maker(atoms, spacer, a, c, xw=1, zw=1):
    """
    This function takes in an array of atom positions one until cell in the slab direction, therefore infinite direction, with ln stacks
    The function outputs an array of atom positions xw and/or zw unit cells wide, therefore making non-inifinite slabs.
    Note:the atom positions are shifted in the positive direction by 10A to account for padding, only if xw or zw !=0
    :param atoms: atom positions for multilayer
    :param spacer: spacing in the unit cell
    :param a: a lattice length
    :param c: c lattice length
    :param xw: number of unit cells wide ribbon will be
    :param zw: number of unit cells wide ribbon will be
    :return: make ribbons. currently edge types not supported
    """

    # create a np array of all permutations in 3 columns e.g. 1 0 1 etc
    base = repeater(permuter(xw, zw), atoms.shape[0]) # stack it up so the permutations are repeated correctly.
                                                      # each permutation is repeated number of atoms times.
                                                    # A permutation refers to the shift in number of unit cells


    combos = zw * xw # the total number of permutations


    # add buffer spacing at the start
    atoms[:,0] = add_buffer(atoms[:,0], spacer, xw) # add spacing in x and z if we are making ribbons
    atoms[:, 2] = add_buffer(atoms[:,2], spacer, zw) #added as quick fix as have sinced changed code


    atoms = stacker(atoms,combos) # repeat the number of atoms combos time so all the permutations are there


    # add atoms positions to the permutations shifted by the unit cell
    unit_cell_shifts = combine(base, array_builder(a, 0, c)) # multiply the unit cell shifts by unit cell sizes

    atoms = adder(atoms, unit_cell_shifts) # add the shifts in cell positions


    # edit current edge so that the ribbons are spatially symmetric
    atoms, distance_chopped, no_deleted = edge_modifier(atoms, xw, zw, 'zz', 2)  # two z position atoms to delete on edge

    return atoms, distance_chopped, no_deleted



def edge_modifier(atoms, xw=1, zw=1, edge_type='zz', atoms_to_delete=2):
    """
    funciton to delete the very edge atoms so that the edges created for the ribbons are spatially symmetric. i.e.
    both finish at the same height. note if xw and zw==1 then nothing changed
        # make cliff edges, and zigzag edges
        # zigzag edges have to be spatially symmetric
    :param atoms: the np array of all atom positions for ribbon created
    :param xw: the width in unit cells of the ribbon in x
    :param zw: the width in unit cells of the ribbon in z
    :return:
    """

    number_atoms_deleted = 0
    distance_chopped = 0

    if (xw-1) or (zw-1):
        if edge_type == 'zz':

            # zigzag - delete the last two z atoms so is symmetric (equivalent to deleting half a unit cell)

            unique, counts = np.unique(atoms[:, 2], return_counts=True) # return unique z values and their counts
            z_values_dict = dict(zip(unique, counts)) # zip into a dictionary for easy acces


            # iterate through deleting biggest elements
            for i in range(atoms_to_delete):

                max_z = unique [-(i+1)] # value of the biggest element
                index = np.argwhere(atoms[:, 2] == max_z) # return indices of elements with largest z value
                atoms = np.delete(atoms, index[:,0], axis = 0) # delete the whole column of these

                number_atoms_deleted += z_values_dict[max_z] # keep track of how many you are deleting!


            distance_chopped = unique[-1] - unique[-3] # how much distance have we chopped off the unit cell


        else:
            print('only zigzag edge types are supported right now') # fix this in the future


    return atoms, distance_chopped, number_atoms_deleted

def find_layers(atoms, layer_thickness = 2.5):

    atoms_layers = dict()

    unique = np.unique(atoms[:,1]) #unique height values

    #this reduces unique to heights which are within distinct layers only
    k=1
    while k:
        try:
            if unique[k] - layer_thickness <= unique[k-1] <= unique[k] + layer_thickness:
                unique = np.delete(unique, k , 0)
            elif k <= len(unique):
                k+=1

        except IndexError as e:
            break

    # scan through and sort into layers
    for j in unique:
        for i in atoms:
            i = i.tolist()
            if ((j-layer_thickness) <= i[1] <= (j+layer_thickness)) and not any((i == x).all() for x in atoms_layers):
                if j in atoms_layers:
                    atoms_layers[j].append(i)
                else:
                    atoms_layers[j] = [i]
    #now atoms_layers should be a dictionary where the key is a height within a layer, and the values are

    return atoms_layers

    


def add_pass_atoms(atoms , atom_type='' , distance_from_edge = 1, type = 'zz', density=1, edge_sensitivity = 1):

    # if density is 1/n then add atoms every nth atom along the periodic direction
    # Then merge with the edge modifier function to find edge atom positions
    # Then write to csv as for intercalant atoms
    #

    #atoms is a numpy array, start by finding the edges according to the type.
    atoms_layers = find_layers(atoms)
    edge_atoms = np.array([])


    if type=='zz':
        #iterate over layers
        for i in atoms_layers.keys():
            atom_layers_np = np.asarray(atoms_layers[i], dtype = np.float32)

            max_width = np.amax(atom_layers_np[:,2]) # maximum edge value
            min_width = np.amin(atom_layers_np[:,2]) # minimum edge value

            for j in atoms_layers[i]: #iterate over atoms in the layer we are looking at
                if min_width - edge_sensitivity <= j[2] <= min_width + edge_sensitivity:
                    edge_atoms = np.append(edge_atoms, [j[0], j[1], j[2] - distance_from_edge], axis = 0)
                elif max_width - edge_sensitivity <= j[2] <= max_width + edge_sensitivity:
                    edge_atoms = np.append(edge_atoms, [j[0], j[1], j[2] + distance_from_edge], axis=0)

    elif type=='ac':
        print('ac not supported yet')
        return False
    else:
        print('invalid type of ribbon passed')
        return False

    edge_atoms = np.reshape(edge_atoms, (-1,3))

    return edge_atoms



def set_stacking_periods(stacking):

    if stacking == 'AA':
        pera, perb, perc = 1, 0, 1
    elif stacking == 'AB':
        pera, perb, perc = 2, 0, 1
    elif stacking == 'AC':
        pera, perb, perc = 2, 0, 2
    elif stacking == 'AD':
        pera, perb, perc = 1, 0, 2
    else:
        print('USE VALID STACKING!')
        pera, perb, perc = 2, 0, 1

    return pera, perb, perc



################## intercalation ##################################


#step 1 -









####################################################################
######################### FUNC DETAILS END #########################
####################################################################




####################################################################
####################### CALL DETAILS START #########################
####################################################################

ln = 4
xw = 1
zw = 3
print('Layer number: '+str(ln)+'. Width of zz ribbon: '+str(zw))

assert type(ln) == int, 'please insert an integer number of layers'
assert ln > 0, 'please insert a positive number of layers'

assert type(xw) == int, 'please insert an integer number of unit cells width'
assert xw > 0, 'please insert a positive number of unit cells width'

assert type(zw) == int, 'please insert an integer number of unit cells width'
assert zw > 0, 'please insert a positive number of unit cells width'

####################################################################
####################### CALL DETAILS END ###########################
####################################################################


####################################################################
################### MATERIAL DETAILS START #########################
####################################################################

#general details
atom_type = "P"
atoms_per_cell = 4
spacer = 10 # Angstrom spacing on each side for sheets or ribbons

stacking = 'AC' # stacking type
print('stacking type = '+str(stacking))
pera, perb, perc = set_stacking_periods(stacking) # periods for different stacking types

# define unit cell ranges. All defined from converged converged bulk pe-d2 vdw
axes = 3 # i.e.  a,b,c --redudant at the moment but may try to use in the future for more general cases
a = 3.342886
b = 10.459090
c = 4.374831
lattice_const = np.array([[a, b, c]])

# thickness of the layers and gaps etc
bii = 0.602462 - 0.397538 # thickness of one layer as fraction of unit cell
thickness = bii*b # thickness in Angstroms
ygap = (0.897538 - 0.602462) * b # top of one layer to bottom of next
yoffset = (0.897538 - 0.397538) * b # bottom of one layer to bottom of next


# intercalation details
intercalant_atom = ''
intercalant_atom_no = 0
intercalant_atom_positions = np.array([])

####################################################################
################### MATERIAL DETAILS END ###########################
####################################################################





####################################################################
####################### CALL FUNCS START ###########################
####################################################################
def edit_cell():
    # build monolayer
    mono = monolayer_builder(spacer, thickness)

    # stack monolayers
    atoms = atom_stacker(mono, ln, atoms_per_cell, axes, lattice_const, yoffset, stacking, pera, perb, perc)

    # if making ribbons, make them here
    atoms, distance_chopped, no_deleted = ribbon_maker(atoms, spacer, a, c, xw, zw)

    # adjust unit cell size
    unit_cell = cell_expander(a, thickness, c, ln, ygap, spacer, xw, zw, distance_chopped)

    combo_rib_dim = xw * zw  # number of permutations of unit cells

    # create intercalation atom cell
    #intercalant_atom_positions = intercalate(atoms, layers, intercalant_atom, stacking='AB', xwidth=1, zwidth=1)

    #add passivated edges
    edge_atoms = add_pass_atoms(atoms, 'H')
    edge_atom_no = len(edge_atoms)


    #IMPROVE YOUR METHOD. DO THIS WITH A DICTIONARY OF ATOM TYPES:NP ARRAYS

    # create poscar string for the given inputs
    poscar_string = build_poscar(atoms, unit_cell, ln, atom_type, atoms_per_cell, combo_rib_dim, no_deleted , 'H', edge_atom_no, edge_atoms) #, intercalant_atom, intercalant_atom_no, intercalant_atom_positions

    # write the poscar to file with filename
    write_to_file = write_poscar(poscar_string, filename)

    return atoms



filename = str(ln)+'_layers_'+str(stacking)+'_stacking_'+'_ribbonwidth_x_'+str(xw)+'_ribbonwidth_z_'+str(zw)+'_zz_edge'
atoms = edit_cell()
####################################################################
######################### CALL FUNCS END ###########################
####################################################################
