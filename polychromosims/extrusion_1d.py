# Stolen from Max, might be added to polychrom at some point

import numpy as np

class leg(object):
    def __init__(self, pos, attrs={"stalled":False, "CTCF":False}):
        """
        A leg has two important attribues: pos (positions) and attrs (a custom list of attributes)
        """
        self.pos = pos
        self.attrs = dict(attrs)

class cohesin(object):
    """
    A cohesin class provides fast access to attributes and positions 
    
    
    cohesin.left is a left leg of cohesin, cohesin.right is a right leg
    cohesin[-1] is also a left leg and cohesin[1] is a right leg         
    
    Also, cohesin.any("myattr") is True if myattr==True in at least one leg
    cohesin.all("myattr") is if myattr=True in both legs
    """
    def __init__(self, leg1, leg2):
        self.left = leg1
        self.right = leg2
   
    def any(self, attr):
        return self.left.attrs[attr] or self.right.attrs[attr]
    
    def all(self, attr):
        return self.left.attrs[attr] and self.right.attrs[attr]    
    
    def __getitem__(self, item):
        if item == -1:
            return self.left
        elif item == 1:
            return self.right 
        else:
            raise ValueError()
        

def unloadProb(cohesin, occupied, args):
    """
    Defines unload probability based on a state of cohesin 
    """
    if occupied[cohesin.left.pos-1] == 2 or occupied[cohesin.right.pos+1] == 2:
        return 1

    if cohesin.any("stalled"):
        # if one side is stalled, we have different unloading probability 
        # Note that here we define stalled cohesins as those stalled not at CTCFs 
        return 1 / args["LIFETIME_STALLED"]
    # otherwise we are just simply unloading 
    return 1 / args["LIFETIME"]    
    


def loadOne(cohesins, occupied, args): 
    """
    A function to load one cohesin 
    """
    while True:
        a = np.random.randint(args["N"])
        if (occupied[a] == 0) and (occupied[a+1] == 0):
            occupied[a] = 1
            occupied[a+1] = 1 
            cohesins.append(cohesin(leg(a), leg(a+1)))
            break


def capture(cohesin, occupied, args):
    """
    We are describing CTCF capture here. 
    This function is specific to this particular project, and 
    users are encouraged to write functions like this 
    
    Note the for-loop over left/right sites below, and using cohesin[side] 
    to get left/right leg. 
    
    Also note how I made ctcfCapture a dict with -1 coding for left side, and 1 for right side 
    and ctcfCapture are dicts as well: keys are locations, and values are probabilities of capture
    """    
    for side in [1, -1]:
        # get probability of capture or otherwise it is 0 
        if np.random.random() < args["ctcfCapture"][side].get(cohesin[side].pos, 0):  
            cohesin[side].attrs["CTCF"] = True  # captured a cohesin at CTCF     
    return cohesin 


def release(cohesin, occupied, args):
    
    """
    AN opposite to capture - releasing cohesins from CTCF 
    """
    
    if not cohesin.any("CTCF"):
        return cohesin  # no CTCF: no release necessary 
        
    # attempting to release either side 
    for side in [-1, 1]: 
        if (np.random.random() < args["ctcfRelease"][side].get(cohesin[side].pos, 0)) and (cohesin[side].attrs["CTCF"]):
            cohesin[side].attrs["CTCF"] = False 
    return cohesin 


def translocate(cohesins, occupied, args):
    """
    This function describes everything that happens with cohesins - 
    loading/unloading them and stalling against each other 
    
    It relies on the functions defined above: unload probability, capture/release. 
    """
    # first we try to unload cohesins and free the matching occupied sites 
    for i in range(len(cohesins)):
        prob = unloadProb(cohesins[i], occupied, args)
        if np.random.random() <= prob:
            occupied[cohesins[i].left.pos] = 0 
            occupied[cohesins[i].right.pos] = 0 
            del cohesins[i]
            loadOne(cohesins, occupied, args)
    
    # then we try to capture and release them by CTCF sites 
    for i in range(len(cohesins)):
        cohesins[i] = capture(cohesins[i], occupied, args)
        cohesins[i] = release(cohesins[i], occupied, args)
    
    # finally we translocate, and mark stalled cohesins because 
    # the unloadProb needs this 
    for i in range(len(cohesins)):
        cohesin = cohesins[i] 
        for leg in [-1,1]: 
            if not cohesin[leg].attrs["CTCF"]: 
                # cohesins that are not at CTCFs and cannot move are labeled as stalled 
                if occupied[cohesin[leg].pos  + leg] != 0:
                    cohesin[leg].attrs["stalled"] = True
                else:
                    cohesin[leg].attrs["stalled"] = False 
                    occupied[cohesin[leg].pos] = 0
                    occupied[cohesin[leg].pos + leg] = 1
                    cohesin[leg].pos += leg        
        cohesins[i] = cohesin

def run_1d(N_mono, N_SMC, frames, CTCFs, lifetime=100, p_capture=0.5,
        p_release=0.02, lifetime_stalled=None, boundaries={1:[],2:[]}):
    """
    Do the 1d simulation

    Parameters
    ----------
    N_mono : int
        number of monomers (i.e. sites in the simulation)
    N_SMC : int
        number of SMCs to throw on
    frames : int
        number of steps to simulate
    CTCFs : array-like
        indices of CTCF sites
    p_capture : float
        probability that a CTCF will catch a passing cohesin
    p_release : float
        probability per step that a CTCF bound cohesin leg will escape
    lifetime : int
        mean lifetime of cohesins, in steps
    lifetime_stalled : int
        mean lifetime of cohesins when stalled against each other. None by
        default, in which case it will be equal to lifetime
    boundaries : dict of lists of int
        set structure of the simulation. boundaries[key] should be a list of
        indices at which to insert key into the 'occupied' array. Use 1 for
        hard boundaries, 2 for strand ends (cohesin unloads upon encounter)

    Returns
    -------
    LEFpositions : np.array, (frames, N_SMC, 2)
        the positions of the LEFs' left and right arms at each time step
    """
    if lifetime_stalled is None:
        lifetime_stalled = lifetime

    ctcfCapture = {ctcf : p_capture for ctcf in CTCFs}
    ctcfRelease = {ctcf : p_release for ctcf in CTCFs}

    args = {}
    args['ctcfRelease'] = {-1 : ctcfRelease, 1 : ctcfRelease}
    args['ctcfCapture'] = {-1 : ctcfCapture, 1 : ctcfCapture}
    args['N'] = N_mono
    args['LIFETIME'] = lifetime
    args['LIFETIME_STALLED'] = lifetime_stalled

    occupied = np.zeros(N_mono)
    for key in boundaries.keys():
        occupied[boundaries[key]] = key
    if occupied[0] == 0:
        occupied[0] = 2
        print('Warning: no boundary specified at 0. Assuming strand end.')
    if occupied[-1] == 0:
        occupied[-1] = 2
        print('Warning: no boundary specified at end. Assuming strand end.')

    cohesins = []
    for _ in range(N_SMC):
        loadOne(cohesins, occupied, args)

    LEFpositions = []
    for _ in range(frames):
        translocate(cohesins, occupied, args)
        LEFpositions.append([(coh.left.pos, coh.right.pos) for coh in cohesins])

    return np.array(LEFpositions)
