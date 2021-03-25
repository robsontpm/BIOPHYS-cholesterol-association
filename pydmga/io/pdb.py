# coding: utf-8
"""
This module is an attempt to grant very basic access to PDB files containing MD trajectories.
Currently it provides one function: readnext(), which returns next line as a tokenized tuple of values,
depending on the type of the PDB record.

Currently supported tokens: "MODEL", "ATOM", "CRYST1", "ENDMDL", "TER"

If token is not supported then function returns tuple with two values. First is token and the second is string "NOT SUPPORTED YET"

This module will be eventually expanded in the future. 
"""
        

def _pdb_convert_number(func, literal):
    '''
    helper function, apply func on literal.strip() and returns its value or 0 if literal.strip() == ''
    func usually is int(), float() etc... 
    '''
    literal = literal.strip()
    return func(literal) if literal else 0


def parserecord(line):
    '''
    given a line of text tries to parse it to PDB record in tuple form
    i.e. (RECORD_NAME, ...). If the line cannot be parset then it returns None
    '''
    while len(line) < 80:
        line = line + ' '
    record_name = line[0:6].strip()
    
    if record_name == "CRYST1":
        #  1 -  6       Record name    "CRYST1"
        #  7 - 15       Real(9.3)      a (Angstroms)
        # 16 - 24       Real(9.3)      b (Angstroms)     
        # 25 - 33       Real(9.3)      c (Angstroms)     
        # 34 - 40       Real(7.2)      alpha (degrees)   
        # 41 - 47       Real(7.2)      beta (degrees)    
        # 48 - 54       Real(7.2)      gamma (degrees)   
        # 56 - 66       LString        Space group       
        # 67 - 70       Integer        Z value             
        a = _pdb_convert_number(float, line[6:15].strip())
        b = _pdb_convert_number(float, line[15:24].strip())
        c = _pdb_convert_number(float, line[24:33].strip())
        alpha = _pdb_convert_number(float, line[33:40].strip())
        beta = _pdb_convert_number(float, line[40:47].strip())
        gamma = _pdb_convert_number(float, line[47:54].strip())
        group = line[55:66]
        z_value = _pdb_convert_number(int, line[66:70].strip())
        return (record_name, a, b, c, alpha, beta, gamma, group, z_value)
            
    if record_name[0:5] == "ORIGX":
        # in record_name[5] we have 1 2 or 3
        #  1 -  6       Record name     "ORIGXn" (n=1, 2, or 3)
        # 11 - 20       Real(10.6)      o[n][1]  
        # 21 - 30       Real(10.6)      o[n][2]  
        # 31 - 40       Real(10.6)      o[n][3]  
        # 46 - 55       Real(10.5)      t[n]             
        return (record_name, "NOT SUPPORTED YET")
    
    if record_name[0:5] == "SCALE":
        # in record_name[5] we have 1 2 or 3
        #  1 -  6       Record name    "SCALEn" (n=1, 2, or 3)
        # 11 - 20       Real(10.6)     s[n][1]                         
        # 21 - 30       Real(10.6)     s[n][2]                       
        # 31 - 40       Real(10.6)     s[n][3]                         
        # 46 - 55       Real(10.5)     u[n]               
        return (record_name, "NOT SUPPORTED YET") 
       
    if record_name[0:5] == "MTRIX":
        # in record_name[5] we have 1 2 or 3
        #  1 -  6       Record name    "MTRIXn" (n=1, 2, or 3)
        #  8 - 10       Integer        Serial number             
        # 11 - 20       Real(10.6)     m[n][1]     
        # 21 - 30       Real(10.6)     m[n][2]     
        # 31 - 40       Real(10.6)     m[n][3]     
        # 46 - 55       Real(10.5)     v[n]        
        # 60            Integer        1 if coordinates for the related molecule are present;
        #                              otherwise, blank.            
        return (record_name, "NOT SUPPORTED YET")
    
    if record_name[0:5] == "TVECT":            
        # in record_name[5] we have space, so it would be stripped away! We can use record_name plainly
        #  1 -  6       Record name    "TVECT "                                    
        #  8 - 10       Integer        Serial number
        # 11 - 20       Real(10.5)     t[1]
        # 21 - 30       Real(10.5)     t[2]
        # 31 - 40       Real(10.5)     t[3]
        # 41 - 70       String         Text comment                        
        return (record_name, "NOT SUPPORTED YET")
    
    if record_name[0:5] == "MODEL":
        # in record_name[5] we have space, so it would be stripped away! We can use record_name plainly
        #  1 -  6       Record name    "MODEL "                                            
        # 11 - 14       Integer        Model serial number
        serial_no = _pdb_convert_number(int, line[10:14]);
        return (record_name, serial_no)
    
    if record_name[0:4] == "ATOM" or record_name == "HETATM":
        # in record_name[4:6] we have spaces, so it would be stripped away! We can use record_name plainly
        #  1 -  6        Record name     "HETATM" or "ATOM  "                                       
        #  7 - 11        Integer         Atom serial number.                
        # 13 - 16        Atom            Atom name                         
        # 17             Character       Alternate location indicator      
        # 18 - 20        Residue name    Residue name                      
        # 22             Character       Chain identifier                  
        # 23 - 26        Integer         Residue sequence number           
        # 27             AChar           Code for insertion of residues    
        # 31 - 38        Real(8.3)       Orthogonal coordinates for X      
        # 39 - 46        Real(8.3)       Orthogonal coordinates for Y      
        # 47 - 54        Real(8.3)       Orthogonal coordinates for Z      
        # 55 - 60        Real(6.2)       Occupancy                         
        # 61 - 66        Real(6.2)       Temperature factor                
        # 73 - 76        LString(4)      Segment identifier, left-justified                    
        # 77 - 78        LString(2)      Element symbol, right-justified  
        # 79 - 80        LString(2)      Charge on the atom              
        serial_no = _pdb_convert_number(int, line[6:11].strip())
        atom_name = line[12:16].strip()
        alternate_location_indicator = line[16]
        residue_name = line[17:20].strip()
        chain = line[21]
        residue_seq_no = _pdb_convert_number(int, line[22:26].strip())
        insertion_code = line[26]
        x = _pdb_convert_number(float, line[30:38].strip())
        y = _pdb_convert_number(float, line[38:46].strip())
        z = _pdb_convert_number(float, line[46:54].strip())
        r = _pdb_convert_number(float, line[54:60].strip())
        temp_factor = _pdb_convert_number(float, line[60:66].strip())
        segment_id = line[72:76].strip()
        element_symbol = line[76:78].strip()
        charge = line[78:80].strip()
        return (record_name,
                serial_no,
                atom_name,
                alternate_location_indicator,
                residue_name,
                chain,
                residue_seq_no,
                insertion_code,
                x, y, z, r,
                temp_factor,
                segment_id,
                element_symbol,
                charge)
    
    if record_name == "ANISOU":
        #  1 -  6        Record name     "ANISOU"                                  
        #  7 - 11        Integer          Atom serial number
        # 13 - 16        Atom             Atom name                  
        # 17             Character        Alternate location indicator                  
        # 18 - 20        Residue name     Residue name               
        # 22             Character        Chain identifier           
        # 23 - 26        Integer          Residue sequence number    
        # 27             AChar           Insertion code             
        # 29 - 35        Integer         u[1][1] 
        # 36 - 42        Integer         u[2][2] 
        # 43 - 49        Integer         u[3][3] 
        # 50 - 56        Integer         u[1][2] 
        # 57 - 63        Integer         u[1][3] 
        # 64 - 70        Integer         u[2][3] 
        # 73 - 76        LString(4)      Segment identifier, left-justified
        # 77 - 78        LString(2)      Element symbol, right-justified
        # 79 - 80        LString(2)      Charge on the atom               
        return (record_name, "NOT SUPPORTED YET")
    
    if record_name[0:4] == "TER":
        # in record_name[3:6] we have spaces, so it would be stripped away! We can use record_name plainly           
        #  1 -  6         Record name       "TER   "                                 
        #  7 - 11         Integer           Serial number
        # 18 - 20         Residue name      Residue name               
        # 22              Character         Chain identifier           
        # 23 - 26         Integer           Residue sequence number    
        # 27              AChar             Insertion code     
        serial_no = _pdb_convert_number(int, line[6:11].strip())
        residue_name = line[17:20].strip()
        chain = line[21]
        residue_seq_no = _pdb_convert_number(int, line[22:26].strip())
        insertion_code = line[26]
        return (record_name, serial_no, residue_name, chain, residue_seq_no, insertion_code)
    
    if record_name == "ENDMDL":
        #  1 -  6         Record name      "ENDMDL"    
        return (record_name,)

    if record_name == "REMARK":
        #  1 -  6         Record name      "REMARK"    (COMMENT)
        return (record_name, line[6:].strip())

    if record_name[0:6] == "TITLE":
        # in record_name[6] we have space, so it would be stripped away! We can use record_name plainly           
        #  1 -  5        Record name      "TITLE"
        tpos = line.find("t=");
        if tpos != -1:
            try:
                time = float(line[tpos+2:].strip())
            except Exception as e:
                str_time = line[tpos+2:].strip().split()[0]            
                time = float(str_time)
        else:
            time = None
        return (record_name, line[6:].strip(), time)        

    return None


def tokenize(line):
    '''
    This is for naming conventions (some our projects used that)
    '''
    return parserecord(line)


def readnext(file, accept_tokens = ["MODEL", "ATOM", "CRYST1", "ENDMDL", "TER", "REMARK", "TITLE"]):
    '''
    it can read normally formatted PDB (that is, the ones compliant with the docs
    (e.g. http://deposit.rcsb.org/adit/docs/pdb_atom_format.html)
    
    :param file: file to read from
    :param accept_tokens: tells which tokens to parse and return, other tokens will be omitted automatically, default = ["MODEL", "ATOM", "CRYST1", "ENDMDL", "TER"]
    :return: tuple representing given token
    
    :exception: StopIteration when EOF approached

    :warning: Currently supported tokens: "MODEL", "ATOM", "CRYST1", "ENDMDL", "TER", "REMARK", "TITLE"

    :warning: If token is not supported then function returns tuple with two values. First is token and the second is string "NOT SUPPORTED YET"
    '''
    while True:
        line = file.readline()        
        if not line:
            raise StopIteration

        record = parserecord(line)
        if record and record[0] in accept_tokens: return record

        
class PDBRecord:
    '''
    TODO: This probably would be better in more OOP way - with PDBRecord - base class and derived classes for various kinds of tokens...

    simple class to represent PDB entry

    it is created on the basic of some tuple, so the simplest way is to pass to it readnext() result
    '''
    def __init__(self, tuple):
        self.tuple = tuple
        self._type = tuple[0]
        # type is allways in tuple[0]
        # and it is stripped string (no redundant spaces!)
        self.locate = {
            'type': 0,
            'a': None,
            'b': None,
            'c': None,            
            'alpha': None,
            'beta': None,
            'gamma': None,
            'space': None,
            'z_value': None,            
            'id': None,
            'atom': None,
            'alt_loc': None,
            'residue': None,
            'chain': None,            
            'residue_seq_no': None,
            'ins_code': None,
            'x': None,
            'y': None,
            'z': None,
            'r': None,
            'temp_fac': None,
            'segment': None,
            'symbol': None,
            'charge': None,
            'text': None, # only for TITLE and REMARK
            'time': None, # only for TITLE
        }
        if self._type == "ATOM" or self._type == "HETATM":
            self.locate['id'] = 1  
            self.locate['atom'] = 2
            self.locate['alt_loc'] = 3
            self.locate['residue'] = 4
            self.locate['chain'] = 5
            self.locate['residue_seq_no'] = 6
            self.locate['ins_code'] = 7
            self.locate['x'] = 8
            self.locate['y'] = 9
            self.locate['z'] = 10
            self.locate['r'] = 11                        
            self.locate['temp_fac'] = 12
            self.locate['segment'] = 13
            self.locate['symbol'] = 14
            self.locate['charge'] = 15          
        if self._type == "CRYST1":
            self.locate['a'] = 1  
            self.locate['b'] = 2
            self.locate['c'] = 3
            self.locate['alpha'] = 4
            self.locate['beta'] = 5
            self.locate['gamma'] = 6
            self.locate['space'] = 7
            self.locate['z_value'] = 8            
        if self._type == "REMARK":
            self.locate['text'] = 1 
        if self._type == "TITLE":
            self.locate['text'] = 1             
            self.locate['time'] = 2                         

    def type(self):
        '''
        returns type of the record, i.e. the first column in PDB file
        '''
        return self._type

    def __getitem__(self, item):
        '''
        gets element of the record by name or None if such item doesn't exist

        available items (depending on the type of the entry):
            
            * 'type'

            * 'a'

            * 'b'

            * 'c'

            * 'alpha'

            * 'beta'

            * 'gamma'

            * 'space'

            * 'z_value'

            * 'id'

            * 'atom'

            * 'alt_loc'

            * 'residue'

            * 'chain'

            * 'residue_seq_no'

            * 'ins_code'

            * 'x'

            * 'y'

            * 'z'

            * 'r'

            * 'temp_fac'

            * 'segment'

            * 'symbol'

            * 'charge'

        '''
        if self.locate[item]:
            return self.tuple[self.locate[item]]
        return None

    def as_coords(self):
        '''
        if this entry has coordinates then it returns tuple (x,y,z), 
        if not then None
        '''
        if self.locate['x']:
            x = self.tuple[self.locate['x']]
            y = self.tuple[self.locate['y']]
            z = self.tuple[self.locate['z']]            
            return (x, y, z)
        return None

    def as_ball(self):
        '''
        if this entry has coordinates and radius then it returns tuple (x,y,z,r), 
        if not then None
        '''
        if self.locate['x'] and self.locate['r']:
            x = self.tuple[self.locate['x']]
            y = self.tuple[self.locate['y']]
            z = self.tuple[self.locate['z']]  
            r = self.tuple[self.locate['r']]
            return (x, y, z, r)        

    def as_particle(self):
        '''
        if this entry has coordinates then it returns tuple (id,x,y,z,r), 
        if not then None
        '''        
        if self.locate['x'] and self.locate['r'] and self.locate['id']:
            id = self.tuple[self.locate['id']] 
            x = self.tuple[self.locate['x']]
            y = self.tuple[self.locate['y']]
            z = self.tuple[self.locate['z']]  
            r = self.tuple[self.locate['r']]
            return (id, x, y, z, r)                 

    def __str__(self):
        '''
        :return: str(this.tuple)
        '''
        return self.tuple.__str__()


def coord_str(pp):
    '''
    For output, returns string representation of 8 digit float with 
    trailing spaces if necessary (format is 8.3 - 8 significant digits 
    with 3 digits after point)
    '''
    tpp = abs(pp)
    spp = "{0}".format(pp)
    if pp >= 0: spp = " " + spp
    if tpp < 10: spp = " " + spp
    if tpp < 100: spp = " " + spp
    while len(spp) < 8: spp = spp + "0"
    return spp[:8]

def radius_str(pp):
    '''
    For output, returns 6 digit float with trailing spaces if necessary 
    (format is 6.2 - 6 significant digits with 2 digits after point)
    '''
    tpp = abs(pp)
    spp = "{0}".format(pp)
    if tpp < 10: spp = " " + spp
    if tpp < 100: spp = " " + spp
    while len(spp) < 6: spp = spp + "0"
    return spp[:6]

def integer_str(pp, digits):
    '''
    For output, returns string representation of the integer exactly digits long, 
    with trailing spaces if necessary.
    '''
    spp = "{0}".format(int(pp)) 
    while len(spp) < digits: spp = " " + spp
    return spp[:digits] 

def hetatm_line(id, x, y, z, r):
    '''
    For output, returns PDB line entry for HETATM (dummy atom)
    '''
    s_id = integer_str(id, 5)
    s_x = coord_str(x)
    s_y = coord_str(y)
    s_z = coord_str(z)
    s_r = radius_str(r)
    return "HETATM{0}                   {1}{2}{3}{4}                   \n".format(s_id, s_x, s_y, s_z, s_r)

def conect_line(atom1_id, atom2_id):
    '''
    For output, returns PDB line entry that forcefuly connects two atoms.
    Usually used to connect HETATMs
    '''    
    s_atom1_id = integer_str(atom1_id, 5)
    s_atom2_id = integer_str(atom2_id, 5)
    return "CONECT{0}{1}               \n".format(s_atom1_id, s_atom2_id)

#
# TODO: kiedys wcielic w życie, ale najpierw poprawic parsowanie PDB, zaby bylo bardziej OOP i latwiejsze w utrzymaniu!
#
# class PDBFrame:
#     def __init__(data, number=None, lines=None):
#         '''
#         number      is the number of the frame.        
#         pdbrecords  are all records read from the file. 
#                     They are of class PDBRecord
#         '''
#         self.lines = lines
#         self.atoms = []        
#         for item in pdbrecords:
#             if item.type() == "ATOM":
#                 self.atoms.append(item)
#             elif item.type() in ["ENDMDL", "TER"]:
#                 self.footer.append(item)
#             elif item["time"] != 
#             else:
#                 self.header.append(item)

#     def make_from_records(self, pdbrecords):

#     def make_from_lines(self, lines):



#     def pdb_str(self):

