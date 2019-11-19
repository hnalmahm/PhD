######################################
###@ Python Script / Crosslinking @###
### Molecular Dynamics Simulations ###
### Using LAMMPS for EPON862-DETDA ###
######################################
###@ Hashim Al Mahmud
###@ MEEM-PhD Student 
###@ Michigan Technological University
###@ Spring 2017 _ updated @ Feb 22, 2017
# email: hnalmahm@mtu.edu
######################################
print 'PYTHON: Here we go!' 
import read
import write

class Bond: pass #.type .atomids = [atom1id, atom2id]
class Angle: pass #.type .atomids = [atom1id, atom2id, atom3id]
class Dihedral: pass #.type .atomids = [atom1id, atom2id, atom3id, atom4id]

m = read.Molecule_File('data.xlinkNC')

a = m.atoms
w = m.masses
b = m.bonds
ang = m.angles
dih = m.dihedrals
b1 = b       # auxiliary dictionaries
b2 = {} 
ang1 = ang
ang2 = {}
dih1 = dih
dih2 = {}

m.extra_lines = "2 extra bond per atom\n"
###########################################################################
################@  step 1: UPDATING BONDS DICTIONARY  @####################
# step 1:A                                                               ## 
# when N (a.type=21) has only one C (a.type=25) new bond N-C (b.type=24) ##
# (1) Create O-H bonds and save it in a new dict b1                      ##
# (2) assign N-H bonds by zero (0) preparing it to delete in another step##
# (3) assign C-O bonds by zero (0) preparing it to delete in another step##
###########################################################################
kk = len(b)
newbonds = 0                                  	 # counter for the action of creating a new bonds
delbonds = 0                                     # counter for the action of deleting a new bonds
N21b = 0                                      	 # counter for the Number of N21 new bonds
for i in b.keys():                        	   	 # locking for N-C New bond
 if b[i].type==24:                   	           # N-C New bond type 24
  if w[a[b[i].atomids[0]].type]==12.011:         # C mass w=12.011
   c25id = b[i].atomids[0]                       # C 25 id
   n21id = b[i].atomids[1]                       # N 21 id
  elif w[a[b[i].atomids[0]].type]==14.007:       # N mass w=14.007
   c25id = b[i].atomids[1]
   n21id = b[i].atomids[0]
  if (a[n21id].type==21):                        # N type 21  in -NC bond
   for j in b.keys():               	           # locking for N-H bond type 17
    if (b[j].type==17):                            # only one of the 2 H45-N21 bonds need to be broken for each new CN bond
     if (b[j].atomids[0]==n21id) or (b[j].atomids[1]==n21id):
      if w[a[b[j].atomids[0]].type]==1.008:       # H mass w=1.008
       h45id = b[j].atomids[0]                    
      elif w[a[b[j].atomids[1]].type]==1.008:    
       h45id = b[j].atomids[1]
      for k in b.keys():                          # O-C bonds type in ether groups
       if b[k].type==1 or b[k].type==9:           # O1-C25 bonds type in ether groups
        if b[k].atomids[0]==c25id or b[k].atomids[1]==c25id:   # C25 same in loop i
         if w[a[b[k].atomids[0]].type]==15.999:    # O mass w=15.999
          o1id = b[k].atomids[0]                    
         elif w[a[b[k].atomids[1]].type]==15.999:  # O mass w=15.999
          o1id = b[k].atomids[1]  
         if (a[o1id].type==23):                     # O type in OH bond
          print 'WARNING!! Repeated OH Bond' 
         else:
          kk = kk + 1                        # increase bonds by 1
          newbond = Bond()                   # create new bond
          newbond.type = 26                  # new bond type (26) for O-H added to the bonds
          newbond.atomids = [o1id,h45id]     # [O,H] Ids
          b1[kk] = newbond                   # add new O-H bond to the auxiliary bonds dictionary b1
          b[j].type = 0                      # lable N-H bond by 0 to exclude it and easily delete it in the next step
          b[k].type = 0                      # lable C-O bond by 0 to exclude it and easily delete it in the next step 
          b1[j].type = 0                     # same action in auxiliary bonds dictionary
          b1[k].type = 0                     # same action in auxiliary bonds dictionary
          newbonds = newbonds + 1
          delbonds = delbonds + 2
          N21b = N21b + 1                    
          a[o1id].type=23                    # O type in OH bond
          a[h45id].type=24                   # H type in OH bond
          a[n21id].type=27                  # N27 to be simply recognized in creating new angles and dihedrals ( note: make sure to re-change it to N21) 
          #print newbonds, '+OH, -CO, -NH' 
##############################################################################
# step 1:B                                                                 ###
# when N (a.type=22) has two C (a.type=25) new bonds N-C (b.type=27)       ###
# (1) Create O-H bonds and save it in a new dict b1                        ###
# (2) assign N-H bonds by zero (0) prepering it to delete in another step  ###
# (3) assign C-O bonds by zero (0) prepering it to delete in another step  ###
##############################################################################
N22b = 0                                         # counter for the Number of N22 new bonds
for i in b.keys():                    	   	 # locking for N-C New bond
 if b[i].type==25:                   	         # N-C New bond type 25
  if w[a[b[i].atomids[0]].type]==12.011:         # C mass w=12.011
   c25id = b[i].atomids[0]                       # C 25 id
   n22id = b[i].atomids[1]                       # N 22/21 id  !!!!!!!!!!!
  elif w[a[b[i].atomids[0]].type]==14.007:       # N mass w=14.007
   c25id = b[i].atomids[1]
   n22id = b[i].atomids[0]
  if (a[n22id].type==22):                        # only N type 22  in -NC bond  (where N22=N26)
   for j in b.keys():               	         # locking for N-H bond type 17
    if b[j].type==17:                            # 2 H-N bonds we need only one to be broken for each new N-C type 27 
     if b[j].atomids[0]==n22id or b[j].atomids[1]==n22id:
      if w[a[b[j].atomids[0]].type]==1.008:       # H mass w=1.008
       h45id = b[j].atomids[0]                    
      elif w[a[b[j].atomids[1]].type]==1.008:    
       h45id = b[j].atomids[1]
      for k in b.keys():                          # O-C bonds type in ether groups
       if b[k].type==1 or b[k].type==9:           # O1-C25 bonds type in ether groups
        if b[k].atomids[0]==c25id or b[k].atomids[1]==c25id:   # C25 same in loop i
         if w[a[b[k].atomids[0]].type]==15.999:    # O mass w=15.999
          o1id = b[k].atomids[0]                    
         elif w[a[b[k].atomids[1]].type]==15.999:    # O mass w=15.999
          o1id = b[k].atomids[1]    
         if (a[o1id].type==23):                     # O type in OH bond
          print 'WARNING!! Repeated OH Bond' 
         else:
          kk = kk + 1                        # increase bonds by 1
          newbond = Bond()                   # create new bond
          newbond.type = 26                  # new bond type (26) for O-H added to the bonds
          newbond.atomids = [o1id,h45id]     # [O,H] Ids
          b1[kk] = newbond                   # add new O-H bond to the auxiliary bonds dictionary b1
          b[j].type = 0                      # lable N-H bond by 0 to exclude it and easily delete it in the next step
          b[k].type = 0                      # lable C-O bond by 0 to exclude it and easily delete it in the next step 
          b1[j].type = 0                     # same action in auxiliary bonds dictionary
          b1[k].type = 0                     # same action in auxiliary bonds dictionary
          newbonds = newbonds + 1
          delbonds = delbonds + 2
          N22b = N22b + 1
#          N21b = N21b - 1
          a[o1id].type=23                    # O type in OH bond
          a[h45id].type=24                   # H type in OH bond
          a[n22id].type=26                  #!!! N26 to be simply recognized in creating new angles and dihedrals (note: make sure to re-change it to N22) 
          #print newbonds, '+OH, -CO, -NH' 
###########################################################################
# step 1:C
# updating bonds dictionaris
# deleting N-H and C-O bonds who are labled by zeros
# rearrange/update bonds dicts (b & b1 using b2)
###########################################################################                
kk = 0
for i in list(b1):
 if b1[i].type == 0:
  del b1[i]
for j in list(b1):
 kk=kk+1
 b2[kk] = b1[j]
b = b2
b1 = b2
m.bonds = b2             
m.nbonds = len(b2)
###########################################################################################################2
################@  step 2: UPDATING ANGLES DICTIONARY  @###################
# step 2:A                                                               ## 
# when N (a.type=21) has only one C (a.type=25) new bond N-C (b.type=24) ##
# (1) Look in the new bonds NC/OH and old bond HN bond type 24           ##
# (2) Create 2 new angles N-C-H  angle type 37                           ##
# (3) Create 1 new angle  N-C-C  angle type 38                           ##
# (4) Create 1 new angles C-N-C  angle type 40                           ##
# (5) Create 1 new angles H-N-C  angle type 41                           ## 
# (6) Create 1 new angles H-O-C  angle type 42                           ##
# (7) in an angle group If Hid in an angle = Hid in OH new bond          ##
#     delete/prep to del its HNC,HNH,OCH,OCC (Type=0)                    ##
###########################################################################
kk = len(ang)
N21a = 0                                        #counter for the the number of N21 new angles
newang21 = 0
delang21 = 0
nch21 = 0
ncc21 = 0
cnc21 = 0                                       # C25-N21-C48
hnc21 = 0
hoc21 = 0
for i in b:
 if b[i].type==24:                               # look for new NC bonds type 24
  if w[a[b[i].atomids[0]].type]==12.011:         # C mass w=12.011
   c25id = b[i].atomids[0]                       # C 25 id
   n21id = b[i].atomids[1]                       # N 21 id
  elif w[a[b[i].atomids[0]].type]==14.007:       # N mass w=14.007
   c25id = b[i].atomids[1]
   n21id = b[i].atomids[0]
  if a[n21id].type==27:                          # check the type of N (if 21)
   N21a = N21a + 1
   for j in b:                                   # look for CH & CC /NC bonds have same C /N in loop i bonds
    if b[j].atomids[0]==c25id or b[j].atomids[1]==c25id:  
     if w[a[b[j].atomids[0]].type]==1.008:              # H46 in CH for new NCH angles
      h46id = b[j].atomids[0]
      kk = kk + 1
      newangle = Angle()                          # create new angles
      newangle.type = 37                          # new angle type (37) for N-C-H added to the angles
      newangle.atomids = [n21id,c25id,h46id]      # [N,C,H] Ids
      ang1[kk] = newangle                         # add new N-C-H angle to the auxiliary angles dictionary ang
      newang21 = newang21 + 1
      nch21 = nch21 + 1
      #print nch21,'new NCH' 
     elif w[a[b[j].atomids[1]].type]==1.008:              # H46 in CH for new NCH angles
      h46id = b[j].atomids[1]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 37                          # new angle type (37) for N-C-H added to the angles
      newangle.atomids = [n21id,c25id,h46id]      # [N,C,H] Ids
      ang1[kk] = newangle                         # add new N-C-H angle to the auxiliary angles dictionary ang
      newang21 = newang21 + 1
      nch21 = nch21 + 1
      #print nch21,'new NCH' 
     elif w[a[b[j].atomids[0]].type]==12.011 and a[b[j].atomids[0]].type != 25: # for new NCC angles
      c13id = b[j].atomids[0]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 38                          # new angle type (38) for N-C-H added to the angles
      newangle.atomids = [n21id,c25id,c13id]      # [N,C,C] Ids
      ang1[kk] = newangle                         # add new N-C-C angle to the auxiliary angles dictionary ang
      newang21 = newang21 + 1
      ncc21 = ncc21 + 1
      #print ncc21,'new NCC' 
     elif w[a[b[j].atomids[1]].type]==12.011 and a[b[j].atomids[1]].type != 25:  # for new NCC angle
      c13id = b[j].atomids[1]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 38                          # new angle type (38) for N-C-C added to the angles
      newangle.atomids = [n21id,c25id,c13id]      # [N,C,C] Ids
      ang1[kk] = newangle                         # add new N-C-C angle to the auxiliary angles dictionary ang
      newang21 = newang21 + 1
      ncc21 = ncc21 + 1
      #print ncc21,'new NCC' 
    elif b[j].atomids[0]==n21id or b[j].atomids[1]==n21id:                     # look for C48 to create new C25-N21-C48 angles
     if w[a[b[j].atomids[0]].type]==12.011 and a[b[j].atomids[0]].type != 25:  # for C48 in NC  bonds
      c48id = b[j].atomids[0]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 40                          # new angle type (40) for C-N-C added to the angles
      newangle.atomids = [c25id,n21id,c48id]      # [C,N,C] Ids
      ang1[kk] = newangle                         # add new C-N-C angle to the auxiliary angles dictionary ang
      newang21 = newang21 + 1
      cnc21 = cnc21 + 1
      #print cnc21,'new CNC' 
     elif w[a[b[j].atomids[1]].type]==12.011 and a[b[j].atomids[1]].type != 25:  # for new CNC angles
      c48id = b[j].atomids[1]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 40                          # new angle type (40) for C-N-C added to the angles
      newangle.atomids = [c25id,n21id,c48id]      # [C,N,C] Ids
      ang1[kk] = newangle                         # add new C-N-C angle to the auxiliary angles dictionary ang
      newang21 = newang21 + 1
      cnc21 = cnc21 + 1
      #print cnc21,'new CNC' 
     elif w[a[b[j].atomids[0]].type]==1.008:      # look for the remaining H45 to create new H45-N21-C25 angles
      h45id = b[j].atomids[0]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 41                          # new angle type (41) for H-N-C added to the angles
      newangle.atomids = [h45id,n21id,c25id]      # [H,N,C] Ids
      ang1[kk] = newangle                         # add new H-N-C angle to the auxiliary angles dictionary ang
      newang21 = newang21 + 1
      hnc21 = hnc21 + 1
      #print hnc21,'new HNC' 
     elif w[a[b[j].atomids[1]].type]==1.008:      # look for the remaining H45 to create new H45-N21-C25 angles
      h45id = b[j].atomids[1]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 41                          # new angle type (41) for H-N-C added to the angles
      newangle.atomids = [h45id,n21id,c25id]      # [H,N,C] Ids
      ang1[kk] = newangle                         # add new H-N-C angle to the auxiliary angles dictionary ang
      newang21 = newang21 + 1
      hnc21 = hnc21 + 1
      #print hnc21,'new HNC' 
 elif b[i].type==26:                             # look for new OH bonds type 26 in loop i
  if w[a[b[i].atomids[0]].type]==15.999:         # O mass w=15.999
   o1id = b[i].atomids[0]                        # O 1 id
   h45id = b[i].atomids[1]                       # H 45 id
  elif w[a[b[i].atomids[0]].type]==1.008:        # H mass w=1.008
   o1id = b[i].atomids[1]
   h45id = b[i].atomids[0]
  for j in b:                                    # look for OC bonds have same O1 in loop i bonds
   if b[j].atomids[0]==o1id or b[j].atomids[1]==o1id:  
    if w[a[b[j].atomids[0]].type]==12.011 and a[b[j].atomids[0]].type != 25:       # for new HOC angles
     c13id = b[j].atomids[0]
     for l in b:
       if b[l].type==24:                         	# look for N21-C25 bonds 
        if w[a[b[l].atomids[0]].type]==12.011:   	# C mass w=12.011
         c25id = b[l].atomids[0]                 	# C 25 id
         n21id = b[l].atomids[1]                    	# N 21 id
        elif w[a[b[l].atomids[0]].type]==14.007:    	# N mass w=14.007
  	 c25id = b[l].atomids[1]
  	 n21id = b[l].atomids[0]
        if a[n21id].type==27:                          # check the type of N (if 21)
          for ll in b:
           if b[ll].atomids[0]==c25id or b[ll].atomids[1]==c25id:
            if b[ll].atomids[0]==c13id or b[ll].atomids[1]==c13id:
             kk = kk + 1
             newangle = Angle()                          # create new angle H45-O1-C13
             newangle.type = 42                          # new angle type (42) for H-O-C added to the angles
             newangle.atomids = [h45id,o1id,c13id]       # [H,O,C] Ids
             ang1[kk] = newangle                         # add new N-C-H angle to the auxiliary angles dictionary ang
    	     newang21 = newang21 + 1
     	     hoc21 = hoc21 + 1
             #print hoc21,'new HOC' 
    elif w[a[b[j].atomids[1]].type]==12.011 and a[b[j].atomids[1]].type != 25:     # for new HOC angles
     c13id = b[j].atomids[1]
     for l in b:
       if b[l].type==24:
        if w[a[b[l].atomids[0]].type]==12.011:         # C mass w=12.011
         c25id = b[l].atomids[0]                       # C 25 id
         n21id = b[l].atomids[1]                       # N 21 id
        elif w[a[b[l].atomids[0]].type]==14.007:       # N mass w=14.007
  	 c25id = b[l].atomids[1]
  	 n21id = b[l].atomids[0]
        if a[n21id].type==27:                          # check the type of N (if 21)
          for ll in b:
           if b[ll].atomids[0]==c25id or b[ll].atomids[1]==c25id:
            if b[ll].atomids[0]==c13id or b[ll].atomids[1]==c13id:
             kk = kk + 1
             newangle = Angle()                          # create new angle
             newangle.type = 42                          # new angle type (42) for H-O-C added to the angles
             newangle.atomids = [h45id,o1id,c13id]       # [H,O,C] Ids
             ang1[kk] = newangle                         # add new N-C-H angle to the auxiliary angles dictionary ang
    	     newang21 = newang21 + 1
     	     hoc21 = hoc21 + 1
             #print hoc21,'new HOC' 
#######################################
## Prepare some old angles to delete ##
#######################################
  for k in list(ang):                            # preper HNH/HNC/OCC/OCH angles to delete (ang[k].type=0)
   if a[ang[k].atomids[1]].type==27:             # only the mid N atom in an angle HNH/HNC   
    if ang[k].atomids[0]==h45id or ang[k].atomids[2]==h45id:         # H id has taken to for OH bond 
     ang[k].type = 0
     ang1[k].type = 0
     delang21 = delang21 + 1                     # HNH/HNC
     #print delang21,'del HNH/HNC' 
   elif a[ang[k].atomids[1]].type==25:           # only the mid C atom in an angle OCC/OCH
    if ang[k].atomids[0]==o1id or ang[k].atomids[2]==o1id: 
      for l in b:
       if b[l].type==24:
        if a[b[l].atomids[0]].type==27 or a[b[l].atomids[1]].type==27:   # only if N type = 21
         if b[l].atomids[0]==ang[k].atomids[1] or b[l].atomids[1]==ang[k].atomids[1]:   # same Cid in loop k and only have bond with N21
          ang[k].type = 0
          ang1[k].type = 0
          delang21 = delang21 + 1                     # one OCC/ two OCH
          #print delang21,'del OCC/OCH' 
################@  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  @###################
# step 2:B      Angles for N22/N26                                       ## 
# when N (a.type=22) has two C (a.type=25) new bond N-C (b.type=27)      ##
# (1) Look in the new bonds NC/OH and old bond HN bond type 27           ##
# (2) Create 2/4 new angles N-C-H  angle type 37                         ##
# (3) Create 1/2 new angle  N-C-C  angle type 38                         ##
# (4) Create 1   new angles C25-N-C25  angle type 39 only for N22  (a)   ##
# (5) Create 1/2 new angles C25-N-C48  angle type 40               (b)   ##
# (6) Create 1/2 new angles H-O-C  angle type 42                         ##
# (7) in an angle group If Hid in an angle = Hid in OH new bond          ##
#     delete/prep to del its 2 HNC, 1 HNH, 4 OCH, 2 OCC (Type=0)         ## 
###########################################################################
N22a = 0                                         # counter for the nuber of N22 new bonds (2 bonds per 1 N22)
newang22 = 0
delang22 = 0
nch22 = 0
ncc22 = 0
cnc22a = 0                                       # C25-N21-C25
cnc22b = 0                                       # C25-N21-C48
hoc22 = 0
checkCNC = 0                                       # checker to prevent repeating C25-N22-C25 angles
for i in b:
 if b[i].type==25:                               # look for new NC bonds type 25
  if w[a[b[i].atomids[0]].type]==12.011:         # C mass w=12.011
   c25id = b[i].atomids[0]                       # C 25 id
   n22id = b[i].atomids[1]                       # N 22 id
  elif w[a[b[i].atomids[0]].type]==14.007:       # N mass w=14.007
   c25id = b[i].atomids[1]
   n22id = b[i].atomids[0]
  if a[n22id].type==26:                          # check the type of N (if 22)
   N22a = N22a + 1                                 
   for j in b:                                   # look for CH & CC /NC bonds have same C /N in loop i bonds
    if b[j].atomids[0]==c25id or b[j].atomids[1]==c25id:  
     if w[a[b[j].atomids[0]].type]==1.008:              # H46 in CH for new NCH angles
      h46id = b[j].atomids[0]
      kk = kk + 1
      newangle = Angle()                          # create new angles
      newangle.type = 37                          # new angle type (37) for N-C-H added to the angles
      newangle.atomids = [n22id,c25id,h46id]      # [N,C,H] Ids
      ang1[kk] = newangle                         # add new N-C-H angle to the auxiliary angles dictionary ang
      newang22 = newang22 + 1
      nch22 = nch22 + 1
      #print nch22,'new N22C25H' 
     elif w[a[b[j].atomids[1]].type]==1.008:              # H46 in CH for new NCH angles
      h46id = b[j].atomids[1]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 37                          # new angle type (37) for N-C-H added to the angles
      newangle.atomids = [n22id,c25id,h46id]      # [N,C,H] Ids
      ang1[kk] = newangle                         # add new N-C-H angle to the auxiliary angles dictionary ang
      newang22 = newang22 + 1
      nch22 = nch22 + 1
      #print nch22,'new N22C25H' 
     elif w[a[b[j].atomids[0]].type]==12.011 and a[b[j].atomids[0]].type != 25: # for new NCC angles
      c13id = b[j].atomids[0]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 38                          # new angle type (38) for N-C-H added to the angles
      newangle.atomids = [n22id,c25id,c13id]      # [N,C,C] Ids
      ang1[kk] = newangle                         # add new N-C-C angle to the auxiliary angles dictionary ang
      newang22 = newang22 + 1
      ncc22 = ncc22 + 1
      #print ncc22,'new N22C25C13' 
     elif w[a[b[j].atomids[1]].type]==12.011 and a[b[j].atomids[1]].type != 25:  # for new NCC angle
      c13id = b[j].atomids[1]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 38                          # new angle type (38) for N-C-C added to the angles
      newangle.atomids = [n22id,c25id,c13id]      # [N,C,C] Ids
      ang1[kk] = newangle                         # add new N-C-C angle to the auxiliary angles dictionary ang
      newang22 = newang22 + 1
      ncc22 = ncc22 + 1
      #print ncc22,'new N22C25C13' 
    elif b[j].type==24:                              # for [C25---24---N22---25---C25] new angles %%%%%%%%%%%%%% 
     if w[a[b[j].atomids[0]].type]==12.011:          # C mass w=12.011
      c25idj = b[j].atomids[0]                       # C 25 id
      n22idj = b[j].atomids[1]                       # N 22 id
     elif w[a[b[j].atomids[0]].type]==14.007:        # N mass w=14.007
      c25idj = b[j].atomids[1]
      n22idj = b[j].atomids[0]
     if a[n22idj].type==26:                          # check the type of N (if 22)
      if (n22idj == n22id) and (c25idj != c25id):
       if checkCNC == n22id:
        print 'Preventing creates a repeated C25-N22-C25 angle' 
       else:
        kk = kk + 1
        newangle = Angle()                          # create new angle
        newangle.type = 39                          # new angle type (39) for N-C-C added to the angles
        newangle.atomids = [c25id,n22id,c25idj]     # [C,N,C] Ids
        ang1[kk] = newangle                         # add new C-N-C angle to the auxiliary angles dictionary ang
        newang22 = newang22 + 1
        cnc22a = cnc22a + 1
        #print cnc22a,'new C25N22C25' 
        checkCNC = n22id
    elif b[j].atomids[0]==n22id or b[j].atomids[1]==n22id:                     # look for C48 to create new C25-N22-C48 angles
     if w[a[b[j].atomids[0]].type]==12.011 and a[b[j].atomids[0]].type != 25:  # for C48 in NC  bonds
      c48id = b[j].atomids[0]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 40                          # new angle type (40) for C-N-C added to the angles
      newangle.atomids = [c25id,n22id,c48id]      # [C,N,C] Ids
      ang1[kk] = newangle                         # add new C-N-C angle to the auxiliary angles dictionary ang
      newang22 = newang22 + 1
      cnc22b = cnc22b + 1
      #print cnc22b,'new C25N22C48' 
     elif w[a[b[j].atomids[1]].type]==12.011 and a[b[j].atomids[1]].type != 25:  # for new CNC angles
      c48id = b[j].atomids[1]
      kk = kk + 1
      newangle = Angle()                          # create new angle
      newangle.type = 40                          # new angle type (40) for C-N-C added to the angles
      newangle.atomids = [c25id,n22id,c48id]      # [C,N,C] Ids
      ang1[kk] = newangle                         # add new C-N-C angle to the auxiliary angles dictionary ang
      newang22 = newang22 + 1
      cnc22b = cnc22b + 1
      #print cnc22b,'new C25N22C48' 
 elif b[i].type==26:                             # look for new OH bonds type 26 in loop i
  if w[a[b[i].atomids[0]].type]==15.999:         # O mass w=15.999
   o1id = b[i].atomids[0]                        # O 1 id
   h45id = b[i].atomids[1]                       # H 45 id
  elif w[a[b[i].atomids[0]].type]==1.008:        # H mass w=1.008
   o1id = b[i].atomids[1]
   h45id = b[i].atomids[0]
  for j in b:                                    # look for OC bonds have same O1 in loop i bonds
   if b[j].atomids[0]==o1id or b[j].atomids[1]==o1id:  
    if w[a[b[j].atomids[0]].type]==12.011 and a[b[j].atomids[0]].type != 25:       # for new HOC angles
     c13id = b[j].atomids[0]
     for l in b:
       if b[l].type==25:                         	# look for N22-C25 bonds 
        if w[a[b[l].atomids[0]].type]==12.011:   	# C mass w=12.011
         c25id = b[l].atomids[0]                 	# C 25 id
         n22id = b[l].atomids[1]                    	# N 22 id
        elif w[a[b[l].atomids[0]].type]==14.007:    	# N mass w=14.007
  	 c25id = b[l].atomids[1]
  	 n22id = b[l].atomids[0]
        if a[n22id].type==26:                          # check the type of N (if 22)
          for ll in b:
           if b[ll].atomids[0]==c25id or b[ll].atomids[1]==c25id:
            if b[ll].atomids[0]==c13id or b[ll].atomids[1]==c13id:
             kk = kk + 1
             newangle = Angle()                          # create new angle H45-O1-C13
             newangle.type = 42                          # new angle type (42) for H-O-C added to the angles
             newangle.atomids = [h45id,o1id,c13id]       # [H,O,C] Ids
             ang1[kk] = newangle                         # add new N-C-H angle to the auxiliary angles dictionary ang
    	     newang22 = newang22 + 1
     	     hoc22 = hoc22 + 1
             #print hoc22,'new HOC13' 
    elif w[a[b[j].atomids[1]].type]==12.011 and a[b[j].atomids[1]].type != 25:     # for new HOC angles
     c13id = b[j].atomids[1]
     for l in b:
       if b[l].type==25:
        if w[a[b[l].atomids[0]].type]==12.011:         # C mass w=12.011
         c25id = b[l].atomids[0]                       # C 25 id
         n22id = b[l].atomids[1]                       # N 22 id
        elif w[a[b[l].atomids[0]].type]==14.007:       # N mass w=14.007
  	 c25id = b[l].atomids[1]
  	 n22id = b[l].atomids[0]
        if a[n22id].type==26:                          # check the type of N (if 22)
          for ll in b:
           if b[ll].atomids[0]==c25id or b[ll].atomids[1]==c25id:
            if b[ll].atomids[0]==c13id or b[ll].atomids[1]==c13id:
             kk = kk + 1
             newangle = Angle()                          # create new angle
             newangle.type = 42                          # new angle type (42) for H-O-C added to the angles
             newangle.atomids = [h45id,o1id,c13id]       # [H,O,C] Ids
             ang1[kk] = newangle                         # add new N-C-H angle to the auxiliary angles dictionary ang
    	     newang22 = newang22 + 1
     	     hoc22 = hoc22 + 1
             #print hoc22,'new HOC13' 
##########################################
## Prepare some old angles to delete    ##
##########################################
for i in b:
 if b[i].type==26:                             # look for new OH bonds type 26 in loop i
  if w[a[b[i].atomids[0]].type]==15.999:         # O mass w=15.999
   o1id = b[i].atomids[0]                        # O 1 id
   h45id = b[i].atomids[1]                       # H 45 id
  elif w[a[b[i].atomids[0]].type]==1.008:        # H mass w=1.008
   o1id = b[i].atomids[1]
   h45id = b[i].atomids[0]
  for k in list(ang):                            # preper OCC/OCH angles to delete (ang[k].type=0)
   if a[ang[k].atomids[1]].type==25:             # only the mid C25 atom in an angle OCC/OCH
    if ang[k].atomids[0]==o1id or ang[k].atomids[2]==o1id: 
      for l in b:
       if b[l].type==25:
        if (a[b[l].atomids[0]].type==26) or (a[b[l].atomids[1]].type==26):   # only if N type = 22
         if b[l].atomids[0]==ang[k].atomids[1] or b[l].atomids[1]==ang[k].atomids[1]: # same C25id in loop k and only have bond with N21
          ang[k].type = 0
          ang1[k].type = 0
          delang22 = delang22 + 1                     # one OCC/CCO & two OCH/HCO
          #print delang22,'del OCC/CCO/OCH/HCO' 
   elif (a[ang[k].atomids[0]].type==25) or (a[ang[k].atomids[2]].type==25) :    # only the edge C25 atom in an angle COC
    if ang[k].atomids[1]==o1id: 
      for l in b:
       if b[l].type==25:                  #C25-N22=C25-N26 bonds
        if a[b[l].atomids[0]].type==26 or a[b[l].atomids[1]].type==26:   # only if N type = 22
         if (b[l].atomids[0]==ang[k].atomids[0]) or (b[l].atomids[0]==ang[k].atomids[2])\
              or (b[l].atomids[1]==ang[k].atomids[0]) or (b[l].atomids[1]==ang[k].atomids[2]): # same C25id in loop k and only have bond with N21
          ang[k].type = 0
          ang1[k].type = 0
          delang22 = delang22 + 1                     # C25-OC/CO-C25
          #print delang22,'del C25-OC/CO-C25' 
for k in list(ang):
 if (a[ang[k].atomids[1]].type==26) and ((w[a[ang[k].atomids[0]].type]==1.008)\
         or (w[a[ang[k].atomids[2]].type]==1.008)):  # prep HNC25/C25NH or HNC48/C48NH to del when the mid atom is N26 i.e. N22
  ang[k].type = 0
  ang1[k].type = 0
  delang22 = delang22 + 1                     # HNC25/C25NH or HNC48/C48NH
  #print delang22,'del HNC48/HNC25' 
################################################################## XXXXXXXXXXXXX
# IMPORTANT WARNING!                                          ####
# THIS STEP IS ESSENTIAL TO DELETE REDUNDUNT(REPEATED) ANGLES ####
################################################################## 
repeatedang = 0
for i in ang1:
 for j in ang1:
  if (i != j) and (ang1[i].type==ang1[j].type) and (ang1[i].atomids[0]==ang1[j].atomids[0])\
     and (ang1[i].atomids[1]==ang1[j].atomids[1]) and (ang1[i].atomids[2]==ang1[j].atomids[2]):
   ang1[j].type = 0 
   repeatedang = repeatedang + 1
   print 'WARNING!!! REPEATED ANGLE IS DELETED' 
  elif (i != j) and (ang1[i].type==ang1[j].type) and (ang1[i].atomids[0]==ang1[j].atomids[2])\
     and (ang1[i].atomids[1]==ang1[j].atomids[1]) and (ang1[i].atomids[2]==ang1[j].atomids[0]):
   ang1[j].type = 0 
   repeatedang = repeatedang + 1
   print 'WARNING!!! REPEATED ANGLE IS DELETED' 
###########################################################################
# updating angles dictionaris
# rearrange/update angles dicts (ang & ang1 using ang2)
########################################################################### 
kk = 0
for i in list(ang1):
 if ang1[i].type == 0:
  del ang1[i]
for j in list(ang1):
 kk = kk + 1
 ang2[kk] = ang1[j]
ang = ang2
ang1 = ang2
m.angles = ang2             
m.nangles = len(ang2)
###########################################################################
# step 3:A
## updating Dihedrals dictionary
## I- create:
## 50    #H45-N44(21)-C13(25)-C13         New     1/b  (1 1 per new NC bond type 24)
## 51    #H45-N44(21)-C13(25)-H46         New     2/b
## 52    #N44(21/22)-C13(25)-C13-C13      New     1/b
## 53    #N44(21/22)-C13(25)-C13-H46      New	  1/b
## 54    #N44(21/22)-C13(25)-C13-O(5/20)  New	  1/b
## 55    #C13(25)-N44(21/22)-C48-C48      New	  2/b
## 56    #C48-N44(21/22)-C13(25)-C13      New	  1/b
## 57    #C48-N44(21/22)-C13(25)-H46      New	  2/b
## 58    #H(7/45)-O(5/20)-C13-C13/(C25)   New	  2/b
## 59    #H(7/45)-O(5/20)-C13-H46         New	  1/b
## II- delete: 						
## 1 H45-N44(21)-C48-C48			  2/b
## 2 C13-C13-C13(25)-O1(20-5-opls)		  1/b
## 3 O1(20-5-opls)-C13(25)-C13-H46		  1/b
## 4 H45(45==>7 opls)-C13(25)-O1(20-5-opls)-C13   2/b
## 5 H45(45==>7 opls)-C13-O1(20-5-opls)-C13(25)   1/b
## 6 C13-C13-O1(20-5-opls)-C13(25)		  1/b
###########################################################################
#print 'UPDATING DIHEDRAL DICTIONARY' 
kk = len(dih)
N21d = 0   #counter for the the number N21 atoms
newdih21 = 0
deldih21 = 0
dih5021 = 0   # type 50  (new dihedrals)
dih5121 = 0   # type 51
dih5221 = 0   # type 52
dih5321 = 0   # type 53
dih5421 = 0   # type 54
dih5521 = 0   # type 55
dih5621 = 0   # type 56
dih5721 = 0   # type 57
dih5821 = 0   # type 58
dih5921 = 0   # type 59
d1x21 = 0     # counters for deleted dihedrals
d2x21 = 0
d3x21 = 0
d4x21 = 0
d5x21 = 0
for i in b:
 if b[i].type==24:                               # look for new NC bonds type 24
  dihflag1 = False
  dihflag2 = False
  if w[a[b[i].atomids[0]].type]==12.011:         # C mass w=12.011
   c25id = b[i].atomids[0]                       # C 25 id
   n21id = b[i].atomids[1]                       # N 21 id
  elif w[a[b[i].atomids[0]].type]==14.007:       # N mass w=14.007
   c25id = b[i].atomids[1]
   n21id = b[i].atomids[0] 
  if a[n21id].type==27:                          # check the type of N (if 21)
   N21d = N21d + 1
   for j in b:               # for dih type 50 / 51
    if (b[j].type==17) and ((b[j].atomids[0]==n21id) or (b[j].atomids[1]==n21id)): # 17 is the type of H45-N
     if (w[a[b[j].atomids[0]].type]==1.008):
      h45id = b[j].atomids[0]
      dihflag1 = True
     else:
      h45id = b[j].atomids[1] 
      dihflag1 = True
     for k in b:
      if ((b[k].atomids[0]==c25id) or (b[k].atomids[1]==c25id)) and ((w[a[b[k].atomids[0]].type]==12.011)\
             and (w[a[b[k].atomids[1]].type]==12.011)):   # C25-C13 bond
       if b[k].atomids[0]==c25id:
        c13id = b[k].atomids[1]
        dihflag2 = True
       else:
        c13id = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                          # create new dih
        newdih.type = 50                              
        newdih.atomids = [h45id,n21id,c25id,c13id]   # [H,N,C,C] Ids
        dih1[kk] = newdih                            # add new H-N-C-C angle to the auxiliary dih dictionary  
        newdih21 = newdih21 + 1
        dih5021 = dih5021 + 1
        dihflag2 = False
        #print dih5021,'new dih5021' 
     for k in b: 
      if ((b[k].atomids[0]==c25id) or (b[k].atomids[1]==c25id)) and ((w[a[b[k].atomids[0]].type]==1.008)\
            or (w[a[b[k].atomids[1]].type]==1.008)):   # C25-H46 bond
       if b[k].atomids[0]==c25id:
        h46id = b[k].atomids[1]
        dihflag2 = True
       else:
        h46id = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                          # create new dih
        newdih.type = 51                              
        newdih.atomids = [h45id,n21id,c25id,h46id]   # [H,N,C,H] Ids
        dih1[kk] = newdih                            # add new H-N-C-C angle to the auxiliary dih dictionary  
        newdih21 = newdih21 + 1
        dih5121 = dih5121 + 1
        dihflag2 = False
        #print dih5121,'new dih5121' 
   dihflag1 = False
   dihflag2 = False
   for j in b:      # for dih type 52 / 53 / 54
    if ((b[j].atomids[0]==c25id) or (b[j].atomids[1]==c25id)) and ((w[a[b[j].atomids[0]].type]==12.011)\
             and (w[a[b[j].atomids[1]].type]==12.011)):   # C25-C13 bond
     if a[b[j].atomids[0]].type ==25:
      c13ida = b[j].atomids[1]
      dihflag1 = True
     else:
      c13ida = b[j].atomids[0] 
      dihflag1 = True
     for k in b:    # for type 52
      if ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and ((w[a[b[k].atomids[0]].type]==12.011)\
         and (w[a[b[k].atomids[1]].type]==12.011)) and ((b[k].atomids[0]!=c25id) and (b[k].atomids[1]!=c25id)):   # C13(a)-C13(b) bond
       if b[k].atomids[0]==c13ida:
        c13idb = b[k].atomids[1]
        dihflag2 = True
       else:
        c13idb = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                          # create new dih
        newdih.type = 52                              
        newdih.atomids = [n21id,c25id,c13ida,c13idb]   # [N,C,C,C] Ids
        dih1[kk] = newdih                            # add new angle to the auxiliary dih dictionary  
        newdih21 = newdih21 + 1
        dih5221 = dih5221 + 1
        dihflag2 = False
        #print dih5221,'new dih5221'   
     for k in b:    # for type 53
      if ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and ((w[a[b[k].atomids[0]].type]==1.008)\
            or (w[a[b[k].atomids[1]].type]==1.008)): # C13(a)-H46 bond
       if b[k].atomids[0]==c13ida:
        h46id = b[k].atomids[1]
        dihflag2 = True
       else:
        h46id = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                          # create new dih
        newdih.type = 53                              
        newdih.atomids = [n21id,c25id,c13ida,h46id]   # [N,C,C,H] Ids
        dih1[kk] = newdih                            # add new angle to the auxiliary dih dictionary  
        newdih21 = newdih21 + 1
        dih5321 = dih5321 + 1
        dihflag2 = False
        #print dih5321,'new dih5321'  
     for k in b:    # for type 54
      if ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and ((w[a[b[k].atomids[0]].type]==15.999)\
            or (w[a[b[k].atomids[1]].type]==15.999)): # C13(a)-O(1) bond
       if b[k].atomids[0]==c13ida:
        o1id = b[k].atomids[1]
        dihflag2 = True
       else:
        o1id = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                          # create new dih
        newdih.type = 54                              
        newdih.atomids = [n21id,c25id,c13ida,o1id]   # [N,C,C,O] Ids
        dih1[kk] = newdih                            # add new angle to the auxiliary dih dictionary  
        newdih21 = newdih21 + 1
        dih5421 = dih5421 + 1
        dihflag2 = False
        #print dih5421,'new dih5421'  
   dihflag1 = False
   dihflag2 = False
   for j in b:       # for dih type 55 / 56 / 57
    if ((b[j].atomids[0]==n21id) or (b[j].atomids[1]==n21id)) and ((w[a[b[j].atomids[0]].type]==12.011)\
             or (w[a[b[j].atomids[1]].type]==12.011)) and (b[j].type != 24):   # N21-C48 bond
     if a[b[j].atomids[0]].type ==27:
      c48ida = b[j].atomids[1]
      dihflag1 = True
     else:
      c48ida = b[j].atomids[0] 
      dihflag1 = True
     for k in b:    # for type 55
      if ((b[k].atomids[0]==c48ida) or (b[k].atomids[1]==c48ida)) and (b[k].atomids[0] != n21id)\
         and (b[k].atomids[1] != n21id):
       if b[k].atomids[0]==c48ida:
        c48idb = b[k].atomids[1]
        dihflag2 = True
       else:
        c48idb = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                          # create new dih
        newdih.type = 55                              
        newdih.atomids = [c25id,n21id,c48ida,c48idb]   # [C,N,C,C] Ids
        dih1[kk] = newdih                            # add new angle to the auxiliary dih dictionary  
        newdih21 = newdih21 + 1
        dih5521 = dih5521 + 1
        dihflag2 = False
        #print dih5521,'new dih5521' 
     for k in b:    # for type 56
      if ((b[k].atomids[0]==c25id) or (b[k].atomids[1]==c25id)) and ((w[a[b[k].atomids[0]].type]==12.011)\
             and (w[a[b[k].atomids[1]].type]==12.011)) and (b[k].type != 24):   # C25-C13 bond
       if b[k].atomids[0]==c25id:
        c13id = b[k].atomids[1]
        dihflag2 = True
       else:
        c13id = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                           # create new dih
        newdih.type = 56                              
        newdih.atomids = [c48ida,n21id,c25id,c13id]   # [C,N,C,C] Ids
        dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
        newdih21 = newdih21 + 1
        dih5621 = dih5621 + 1
        dihflag2 = False
        #print dih5621,'new dih5621' 
     for k in b:    # for type 57
      if ((b[k].atomids[0]==c25id) or (b[k].atomids[1]==c25id)) and ((w[a[b[k].atomids[0]].type]==1.008)\
             or (w[a[b[k].atomids[1]].type]==1.008)) and (b[k].type != 24):   # C25-H46 bond
       if b[k].atomids[0]==c25id:
        h46id = b[k].atomids[1]
        dihflag2 = True
       else:
        h46id = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                           # create new dih
        newdih.type = 57                              
        newdih.atomids = [c48ida,n21id,c25id,h46id]   # [C,N,C,H] Ids
        dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
        newdih21 = newdih21 + 1
        dih5721 = dih5721 + 1
        dihflag2 = False
        #print dih5721,'new dih5721'   
   dihflag1 = False ; dihflag2 = False ; dihflag3 = False ; dihflag4 = False
   for j in b:       # for dih type 58 / 59 
    if ((b[j].atomids[0]==c25id) or (b[j].atomids[1]==c25id)) and ((w[a[b[j].atomids[0]].type]==12.011)\
             and (w[a[b[j].atomids[1]].type]==12.011)) and (b[j].type != 24):   # C25-C13a bond 
     if b[j].atomids[0]==c25id:
      c13ida = b[j].atomids[1]
      dihflag1 = True
     else:
      c13ida = b[j].atomids[0] 
      dihflag1 = True
     for k in b:     
      if ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and ((w[a[b[k].atomids[0]].type]==15.999)\
             or (w[a[b[k].atomids[1]].type]==15.999)):   # C13a-O bond
       if b[k].atomids[0]==c13ida:
        o1id = b[k].atomids[1]
        dihflag2 = True
       else:
        o1id = b[k].atomids[0] 
        dihflag2 = True
     for k in b:    
      if (b[k].type==26) and ((b[k].atomids[0]==o1id) or (b[k].atomids[1]==o1id)):   # O-H bond
       if b[k].atomids[0]==o1id:
        h45id = b[k].atomids[1]
        dihflag3 = True
       else:
        h45id = b[k].atomids[0]
        dihflag3 = True 
       if (dihflag1 == True) and (dihflag2 == True) and (dihflag3 == True):
        kk = kk + 1
        newdih = Dihedral()                           # create new dih
        newdih.type = 58                              
        newdih.atomids = [h45id,o1id,c13ida,c25id]   # [H,O,C,C] Ids
        dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
        newdih21 = newdih21 + 1
        dih5821 = dih5821 + 1
        #print dih5821,'new dih5821'  
        for k in b: 
# for type 58                                                                  
         if ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and\
            ((w[a[b[k].atomids[0]].type]==12.011) and (w[a[b[k].atomids[1]].type]==12.011))\
            and ((b[k].atomids[0]!=c25id) and (b[k].atomids[1]!=c25id)): # C13a-C13b bond 
          if b[k].atomids[0]==c13ida:
           c13idb = b[k].atomids[1]
           dihflag4 = True
          else:
           c13idb = b[k].atomids[0]
           dihflag4 = True 
          if (dihflag1 == True) and (dihflag2 == True) and (dihflag3 == True) and (dihflag4 == True):
           kk = kk + 1
           newdih = Dihedral()                           # create new dih
           newdih.type = 58                              
           newdih.atomids = [h45id,o1id,c13ida,c13idb]   # [H,O,C,C] Ids
           dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
           newdih21 = newdih21 + 1
           dih5821 = dih5821 + 1
           dihflag4 = False
           #print dih5821,'new dih5821' 
           dihflag4 = False
# for type 59
         elif ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and\
            ((w[a[b[k].atomids[0]].type]==1.008) or (w[a[b[k].atomids[1]].type]==1.008)): # C13a-H46 bond 
          if b[k].atomids[0]==c13ida:
           h46id = b[k].atomids[1]
           dihflag4 = True
          else:
           h46id = b[k].atomids[0]
           dihflag4 = True 
          if (dihflag1 == True) and (dihflag2 == True) and (dihflag3 == True) and (dihflag4 == True):
           kk = kk + 1
           newdih = Dihedral()                           # create new dih
           newdih.type = 59                              
           newdih.atomids = [h45id,o1id,c13ida,h46id]   # [H,O,C,H] Ids
           dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
           newdih21 = newdih21 + 1
           dih5921 = dih5921 + 1
           dihflag4 = False
           #print dih5921,'new dih5921' 
   dihflag1 = False ; dihflag2 = False ; dihflag3 = False ; dihflag4 = False 
#############################################
## Prepare some old dihedrals to delete    ##
#############################################
## d1 H45-N44(21)-C48-C48			  2/b
## d2 C13-C13-C13(25)-O1(20-5-opls)		  1/b
## d2 O1(20-5-opls)-C13(25)-C13-H46		  1/b
## d3 H46-C13(25)-O1(20-5-opls)-C13               2/b
## d4 H46-C13-O1(20-5-opls)-C13(25)               1/b
## d5 C13-C13-O1(20-5-opls)-C13(25)		  1/b
for i in b:
 if b[i].type==24:                               # look for new NC bonds type 24
  if w[a[b[i].atomids[0]].type]==12.011:         # C mass w=12.011
   c25id = b[i].atomids[0]                       # C 25 id
   n21id = b[i].atomids[1]                       # N 21 id
  elif w[a[b[i].atomids[0]].type]==14.007:       # N mass w=14.007
   c25id = b[i].atomids[1]
   n21id = b[i].atomids[0] 
  if a[n21id].type==27:                          # check the type of N (if 21)
   for j in b:
    if b[j].type==26:                               # look for new OH bonds type 26
     if w[a[b[j].atomids[0]].type]==15.999:         # O mass w=15.999
      o1id = b[j].atomids[0]                        # O id
      h45id = b[j].atomids[1]                       # H id
     elif w[a[b[j].atomids[0]].type]==1.008:        # H mass w=1.008
      o1id = b[j].atomids[1]
      h45id = b[j].atomids[0] 
     for k in dih:
      if ((dih[k].atomids[0]==h45id) and (dih[k].atomids[1]==n21id)) or\
          ((dih[k].atomids[2]==n21id) and (dih[k].atomids[3]==h45id)):  # 1- HNC48C48 or C48C48NH 
       dih[k].type = 0
       dih1[k].type = 0
       d1x21 = d1x21 + 1                     # HNC48C48/C48C48NH
       deldih21 = deldih21 + 1
       #print d1x21,'deldih HNCC' 
      elif ((dih[k].atomids[0]==o1id) and (dih[k].atomids[1]==c25id)) or\
          ((dih[k].atomids[2]==c25id) and (dih[k].atomids[3]==o1id)):  # OCCC or OCCH
       dih[k].type = 0
       dih1[k].type = 0
       d2x21 = d2x21 + 1  
       deldih21 = deldih21 + 1                  
       #print d2x21,'deldih OCCC/OCCH' 
      elif ((dih[k].atomids[1]==o1id) and (dih[k].atomids[2]==c25id)) or\
          ((dih[k].atomids[1]==c25id) and (dih[k].atomids[2]==o1id)):  # HC25OC/COC25H
       dih[k].type = 0
       dih1[k].type = 0
       d3x21 = d3x21 + 1   
       deldih21 = deldih21 + 1                 
       #print d3x21,'deldih HC25OC' 
      elif ((dih[k].atomids[0]==c25id) and (dih[k].atomids[1]==o1id) and (w[a[dih[k].atomids[3]].type]==1.008)) or\
          ((dih[k].atomids[2]==o1id) and (dih[k].atomids[3]==c25id) and (w[a[dih[k].atomids[0]].type]==1.008)):  # HCOC25
       dih[k].type = 0
       dih1[k].type = 0
       d4x21 = d4x21 + 1   
       deldih21 = deldih21 + 1                 
       #print d4x21,'deldih HCOC25' 
      elif ((dih[k].atomids[0]==c25id) and (dih[k].atomids[1]==o1id) and (w[a[dih[k].atomids[3]].type]==12.011)) or\
           ((dih[k].atomids[3]==c25id) and (dih[k].atomids[2]==o1id) and (w[a[dih[k].atomids[0]].type]==12.011)) or\
           ((dih[k].atomids[3]==o1id) and (dih[k].atomids[2]==c25id) and (w[a[dih[k].atomids[0]].type]==12.011)) or\
           ((dih[k].atomids[0]==o1id) and (dih[k].atomids[1]==c25id) and (w[a[dih[k].atomids[3]].type]==12.011)):  # HCOC25/C25OCH/HCC25O/OC25CH
       dih[k].type = 0
       dih1[k].type = 0
       d5x21 = d5x21 + 1    
       deldih21 = deldih21 + 1                
       #print d5x21,'deldih CCOC25' 
###########################################################################
# step 3:B
## updating Dihedrals dictionary
## I- create:
## 50    
## 51   
## 52    #N44(21/22)-C13(25)-C13-C13      New   1/b  (1 per new NC bond type 27)
## 53    #N44(21/22)-C13(25)-C13-H46      New	1/b
## 54    #N44(21/22)-C13(25)-C13-O(5/20)  New	1/b
## 55    #C13(25)-N44(21/22)-C48-C48      New	2/b
## 56    #C48-N44(21/22)-C13(25)-C13      New	1/b
## 57    #C48-N44(21/22)-C13(25)-H46      New	2/b
## 58    #H(7/45)-O(5/20)-C13-C13/(C25)   New	2/b
## 59    #H(7/45)-O(5/20)-C13-H46         New	1/b
## 60    #C13(25)-N44(22)-C13(25)-C13     New	1/b
## 61    #C13(25)-N44(22)-C13(25)-H46     New	2/b
## II- delete: 						
## H45-N44(21)-C48-C48				2/b
## C13-C13-C13(25)-O1(20-5-opls)		1/b
## O1(20-5-opls)-C13(25)-C13-H46		1/b
## H45(45==>7 opls)-C13(25)-O1(20-5-opls)-C13   2/b
## H45(45==>7 opls)-C13-O1(20-5-opls)-C13(25)   1/b
## C13-C13-O1(20-5-opls)-C13(25)		1/b
###########################################################################
N22d = 0   #counter for the the number N22 atoms
checkN22 = 0    # to prevent repeated dih for the same N22
newdih22 = 0
deldih22 = 0
dih5222 = 0   # type 52 (new dihedrals)
dih5322 = 0   # type 53
dih5422 = 0   # type 54
dih5522 = 0   # type 55
dih5622 = 0   # type 56
dih5722 = 0   # type 57
dih5822 = 0   # type 58
dih5922 = 0   # type 59
dih6022 = 0   # type 60
dih6122 = 0   # type 61
d1x22 = 0     # counters for deleted dihedrals
d2x22 = 0
d3x22 = 0
d4x22 = 0
d5x22 = 0
for i in b:
 if b[i].type==25:                               # look for new NC bonds type 25
  dihflag1 = False
  dihflag2 = False
  if w[a[b[i].atomids[0]].type]==12.011:         # C mass w=12.011
   c25id = b[i].atomids[0]                       # C 25 id
   n22id = b[i].atomids[1]                       # N 22 id
  elif w[a[b[i].atomids[0]].type]==14.007:       # N mass w=14.007
   c25id = b[i].atomids[1]
   n22id = b[i].atomids[0] 
  if a[n22id].type==26:                          # check the type of N (if 22)
   N22d = N22d + 1
   for j in b:      # for dih type 52 / 53 / 54
    if ((b[j].atomids[0]==c25id) or (b[j].atomids[1]==c25id)) and ((w[a[b[j].atomids[0]].type]==12.011)\
             and (w[a[b[j].atomids[1]].type]==12.011)):   # C25-C13 bond
     if a[b[j].atomids[0]].type ==25:
      c13ida = b[j].atomids[1]
      dihflag1 = True
     else:
      c13ida = b[j].atomids[0] 
      dihflag1 = True
     for k in b:    # for type 52
      if ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and ((w[a[b[k].atomids[0]].type]==12.011)\
         and (w[a[b[k].atomids[1]].type]==12.011)) and ((b[k].atomids[0]!=c25id) and (b[k].atomids[1]!=c25id)):   # C13(a)-C13(b) bond
       if b[k].atomids[0]==c13ida:
        c13idb = b[k].atomids[1]
        dihflag2 = True
       else:
        c13idb = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                          # create new dih
        newdih.type = 52                              
        newdih.atomids = [n22id,c25id,c13ida,c13idb]   # [N,C,C,C] Ids
        dih1[kk] = newdih                            # add new angle to the auxiliary dih dictionary  
        newdih22 = newdih22 + 1
        dih5222 = dih5222 + 1
        dihflag2 = False
        #print dih5222,'new dih5222'  
     for k in b:    # for type 53
      if ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and ((w[a[b[k].atomids[0]].type]==1.008)\
            or (w[a[b[k].atomids[1]].type]==1.008)): # C13(a)-H46 bond
       if b[k].atomids[0]==c13ida:
        h46id = b[k].atomids[1]
        dihflag2 = True
       else:
        h46id = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                          # create new dih
        newdih.type = 53                              
        newdih.atomids = [n22id,c25id,c13ida,h46id]   # [N,C,C,H] Ids
        dih1[kk] = newdih                            # add new angle to the auxiliary dih dictionary  
        newdih22 = newdih22 + 1
        dih5322 = dih5322 + 1
        dihflag2 = False
        #print dih5322,'new dih5322'  
     for k in b:    # for type 54
      if ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and ((w[a[b[k].atomids[0]].type]==15.999)\
            or (w[a[b[k].atomids[1]].type]==15.999)): # C13(a)-O(1) bond
       if b[k].atomids[0]==c13ida:
        o1id = b[k].atomids[1]
        dihflag2 = True
       else:
        o1id = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                          # create new dih
        newdih.type = 54                              
        newdih.atomids = [n22id,c25id,c13ida,o1id]   # [N,C,C,O] Ids
        dih1[kk] = newdih                            # add new angle to the auxiliary dih dictionary  
        newdih22 = newdih22 + 1
        dih5422 = dih5422 + 1
        dihflag2 = False
        #print dih5422,'new dih5422'   
   dihflag1 = False
   dihflag2 = False  
   for j in b:       # for dih type 55 / 56 / 57
    if ((b[j].atomids[0]==n22id) or (b[j].atomids[1]==n22id)) and ((w[a[b[j].atomids[0]].type]==12.011)\
             or (w[a[b[j].atomids[1]].type]==12.011)) and ((b[j].type != 24) and (b[j].type != 25)):   # N21-C48 bond
     if a[b[j].atomids[0]].type ==26:
      c48ida = b[j].atomids[1]
      dihflag1 = True
     else:
      c48ida = b[j].atomids[0] 
      dihflag1 = True
     for k in b:    # for type 55
      if ((b[k].atomids[0]==c48ida) or (b[k].atomids[1]==c48ida)) and (b[k].atomids[0] != n22id)\
         and (b[k].atomids[1] != n22id):
       if b[k].atomids[0]==c48ida:
        c48idb = b[k].atomids[1]
        dihflag2 = True
       else:
        c48idb = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                          # create new dih
        newdih.type = 55                              
        newdih.atomids = [c25id,n22id,c48ida,c48idb]   # [C,N,C,C] Ids
        dih1[kk] = newdih                            # add new angle to the auxiliary dih dictionary  
        newdih22 = newdih22 + 1
        dih5522 = dih5522 + 1
        dihflag2 = False
        #print dih5522,'new dih5522' 
     for k in b:    # for type 56
      if ((b[k].atomids[0]==c25id) or (b[k].atomids[1]==c25id)) and ((w[a[b[k].atomids[0]].type]==12.011)\
             and (w[a[b[k].atomids[1]].type]==12.011)) and ((b[k].type != 24) and (b[k].type != 25)):   # C25-C13 bond
       if b[k].atomids[0]==c25id:
        c13id = b[k].atomids[1]
        dihflag2 = True
       else:
        c13id = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                           # create new dih
        newdih.type = 56                              
        newdih.atomids = [c48ida,n22id,c25id,c13id]   # [C,N,C,C] Ids
        dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
        newdih22 = newdih22 + 1
        dih5622 = dih5622 + 1
        dihflag2 = False
        #print dih5622,'new dih5622' 
     for k in b:    # for type 57
      if ((b[k].atomids[0]==c25id) or (b[k].atomids[1]==c25id)) and ((w[a[b[k].atomids[0]].type]==1.008)\
             or (w[a[b[k].atomids[1]].type]==1.008)) and ((b[k].type != 24) and (b[k].type != 25)):   # C25-H46 bond
       if b[k].atomids[0]==c25id:
        h46id = b[k].atomids[1]
        dihflag2 = True
       else:
        h46id = b[k].atomids[0]
        dihflag2 = True
       if (dihflag1 == True) and (dihflag2 == True):
        kk = kk + 1
        newdih = Dihedral()                           # create new dih
        newdih.type = 57                              
        newdih.atomids = [c48ida,n22id,c25id,h46id]   # [C,N,C,H] Ids
        dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
        newdih22 = newdih22 + 1
        dih5722 = dih5722 + 1
        dihflag2 = False
        #print dih5722,'new dih5722'  
   dihflag1 = False ; dihflag2 = False ; dihflag3 = False ; dihflag4 = False
   for j in b:       # for dih type 58 / 59 
    if ((b[j].atomids[0]==c25id) or (b[j].atomids[1]==c25id)) and ((w[a[b[j].atomids[0]].type]==12.011)\
             and (w[a[b[j].atomids[1]].type]==12.011)) and ((b[j].type != 24) and (b[j].type != 25)):   # C25-C13a bond 
     if b[j].atomids[0]==c25id:
      c13ida = b[j].atomids[1]
      dihflag1 = True
     else:
      c13ida = b[j].atomids[0] 
      dihflag1 = True
     for k in b:     
      if ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and ((w[a[b[k].atomids[0]].type]==15.999)\
             or (w[a[b[k].atomids[1]].type]==15.999)):   # C13a-O bond
       if b[k].atomids[0]==c13ida:
        o1id = b[k].atomids[1]
        dihflag2 = True
       else:
        o1id = b[k].atomids[0] 
        dihflag2 = True
     for k in b:    # for type 58
      if (b[k].type==26) and ((b[k].atomids[0]==o1id) or (b[k].atomids[1]==o1id)):   # O-H bond
       if b[k].atomids[0]==o1id:
        h45id = b[k].atomids[1]
        dihflag3 = True
       else:
        h45id = b[k].atomids[0]
        dihflag3 = True 
       if (dihflag1 == True) and (dihflag2 == True) and (dihflag3 == True):
        kk = kk + 1
        newdih = Dihedral()                           # create new dih
        newdih.type = 58                              
        newdih.atomids = [h45id,o1id,c13ida,c25id]   # [H,O,C,C] Ids
        dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
        newdih22 = newdih22 + 1
        dih5822 = dih5822 + 1
        #print dih5822,'new dih5822' 
        for k in b:                                                                   
         if ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and\
            ((w[a[b[k].atomids[0]].type]==12.011) and (w[a[b[k].atomids[1]].type]==12.011))\
            and ((b[k].atomids[0]!=c25id) and (b[k].atomids[1]!=c25id)): # C13a-C13b bond 
          if b[k].atomids[0]==c13ida:
           c13idb = b[k].atomids[1]
           dihflag4 = True
          else:
           c13idb = b[k].atomids[0]
           dihflag4 = True 
          if (dihflag1 == True) and (dihflag2 == True) and (dihflag3 == True) and (dihflag4 == True):
           kk = kk + 1
           newdih = Dihedral()                           # create new dih
           newdih.type = 58                              
           newdih.atomids = [h45id,o1id,c13ida,c13idb]   # [H,O,C,C] Ids
           dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
           newdih22 = newdih22 + 1
           dih5822 = dih5822 + 1
           dihflag4 = False
           #print dih5822,'new dih5822' 
           dihflag4 = False
# for type 59
         elif ((b[k].atomids[0]==c13ida) or (b[k].atomids[1]==c13ida)) and\
            ((w[a[b[k].atomids[0]].type]==1.008) or (w[a[b[k].atomids[1]].type]==1.008)): # C13a-H46 bond 
          if b[k].atomids[0]==c13ida:
           h46id = b[k].atomids[1]
           dihflag4 = True
          else:
           h46id = b[k].atomids[0]
           dihflag4 = True 
          if (dihflag1 == True) and (dihflag2 == True) and (dihflag3 == True) and (dihflag4 == True):
           kk = kk + 1
           newdih = Dihedral()                           # create new dih
           newdih.type = 59                              
           newdih.atomids = [h45id,o1id,c13ida,h46id]   # [H,O,C,H] Ids
           dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
           newdih22 = newdih22 + 1
           dih5922 = dih5922 + 1
           dihflag4 = False
           #print dih5922,'new dih5922' 
   dihflag1 = False ; dihflag2 = False ; dihflag3 = False ; dihflag4 = False 
   for j in b:       # for dih type 60 / 61 
    if b[j].type==24: 
     if w[a[b[j].atomids[0]].type]==12.011:          # C mass w=12.011
      c25id2 = b[j].atomids[0]                       # C 25 id ==>C25--24--N22--25--C25  %%%%%%%%%%%%%%
      n22id2= b[j].atomids[1]                        # N 22 id
     elif w[a[b[j].atomids[0]].type]==14.007:        # N mass w=14.007
      c25id2 = b[j].atomids[1]
      n22id2 = b[j].atomids[0] 
     if (n22id2 == n22id) and (c25id2 != c25id):
      if checkN22 == n22id:   # to prevent repeated dih for the same N22
       print '##################  Skip repeated Dih' 
      else:
       dihflag1 = True
       for k in b:
        if ((b[k].atomids[0]==c25id) or (b[k].atomids[1]==c25id)) and ((w[a[b[k].atomids[0]].type]==12.011)\
             and (w[a[b[k].atomids[1]].type]==12.011)) and ((b[k].type != 24) and (b[k].type != 25)):   # C25-C13 Right bond 
         if b[k].atomids[0]==c25id:
          c13id = b[k].atomids[1]
          dihflag2 = True
         else:
          c13id = b[k].atomids[0]
          dihflag2 = True 
         if (dihflag1 == True) and (dihflag2 == True):
          kk = kk + 1
          newdih = Dihedral()                           # create new dih
          newdih.type = 60                              
          newdih.atomids = [c25id2,n22id,c25id,c13id]   # [C,N,C,C] Ids  Right
          dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
          newdih22 = newdih22 + 1
          dih6022 = dih6022 + 1
          checkN22 = n22id           # to prevent repeated dih for the same N22
          dihflag2 = False
          #print dih6022,'new dih6022' 
        elif ((b[k].atomids[0]==c25id2) or (b[k].atomids[1]==c25id2)) and ((w[a[b[k].atomids[0]].type]==12.011)\
             and (w[a[b[k].atomids[1]].type]==12.011)) and ((b[k].type != 24) and (b[k].type != 25)):   # C25-C13 left bond 
         if b[k].atomids[0]==c25id2:
          c13id2 = b[k].atomids[1]
          dihflag2 = True
         else:
          c13id2 = b[k].atomids[0]
          dihflag2 = True 
         if (dihflag1 == True) and (dihflag2 == True):
          kk = kk + 1
          newdih = Dihedral()                           # create new dih
          newdih.type = 60                              
          newdih.atomids = [c25id,n22id,c25id2,c13id2]   # [C,N,C,C] Ids left
          dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
          newdih22 = newdih22 + 1
          dih6022 = dih6022 + 1
          checkN22 = n22id           # to prevent repeated dih for the same N22
          dihflag2 = False
          #print dih6022,'new dih6022'      
        elif ((b[k].atomids[0]==c25id) or (b[k].atomids[1]==c25id)) and ((w[a[b[k].atomids[0]].type]==1.008)\
             or (w[a[b[k].atomids[1]].type]==1.008)):   # C25-H46 Right bond  
         if b[k].atomids[0]==c25id:
          h46id = b[k].atomids[1]
          dihflag2 = True
         else:
          h46id = b[k].atomids[0]
          dihflag2 = True 
         if (dihflag1 == True) and (dihflag2 == True):
          kk = kk + 1
          newdih = Dihedral()                           # create new dih
          newdih.type = 61                              
          newdih.atomids = [c25id2,n22id,c25id,h46id]   # [C,N,C,H] Ids  Right
          dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
          newdih22 = newdih22 + 1
          dih6122 = dih6122 + 1
          checkN22 = n22id           # to prevent repeated dih for the same N22
          dihflag2 = False 
          #print dih6122,'new dih6122' 
        elif ((b[k].atomids[0]==c25id2) or (b[k].atomids[1]==c25id2)) and ((w[a[b[k].atomids[0]].type]==1.008)\
             or (w[a[b[k].atomids[1]].type]==1.008)):   # C25-H46 left bond 
         if b[k].atomids[0]==c25id2:
          h46id2 = b[k].atomids[1]
          dihflag2 = True
         else:
          h46id2 = b[k].atomids[0]
          dihflag2 = True 
         if (dihflag1 == True) and (dihflag2 == True):
          kk = kk + 1
          newdih = Dihedral()                           # create new dih
          newdih.type = 61                              
          newdih.atomids = [c25id,n22id,c25id2,h46id2]   # [C,N,C,C] Ids left
          dih1[kk] = newdih                             # add new angle to the auxiliary dih dictionary  
          newdih22 = newdih22 + 1
          dih6122 = dih6122 + 1
          checkN22 = n22id           # to prevent repeated dih for the same N22
          dihflag2 = False
          #print dih6122,'new dih6122'            
#############################################
## Prepare some old dihedrals to delete    ##
#############################################
## d1 H45-N44(22)-C48-C48			  2/b
## d2 C13-C13-C13(25)-O1(20-5-opls)		  1/b
## d2 O1(20-5-opls)-C13(25)-C13-H46		  1/b
## d3 H46-C13(25)-O1(20-5-opls)-C13               2/b
## d4 H46-C13-O1(20-5-opls)-C13(25)               1/b
## d5 C13-C13-O1(20-5-opls)-C13(25)		  1/b
#############################################
## Prepare some old dihedrals to delete    ##
#############################################
for i in b:
 if b[i].type==25:                               # look for new NC bonds type 25
  if w[a[b[i].atomids[0]].type]==12.011:         # C mass w=12.011
   c25id = b[i].atomids[0]                       # C 25 id
   n22id = b[i].atomids[1]                       # N 22 id
  elif w[a[b[i].atomids[0]].type]==14.007:       # N mass w=14.007
   c25id = b[i].atomids[1]
   n22id = b[i].atomids[0] 
  if a[n22id].type==26:                          # check the type of N (if 22)
   for j in b:
    if b[j].type==26:                               # look for new OH bonds type 26
     if w[a[b[j].atomids[0]].type]==15.999:         # O mass w=15.999
      o1id = b[j].atomids[0]                        # O id
      h45id = b[j].atomids[1]                       # H id
     elif w[a[b[j].atomids[0]].type]==1.008:        # H mass w=1.008
      o1id = b[j].atomids[1]
      h45id = b[j].atomids[0] 
     for k in dih:
      if ((dih[k].atomids[0]==h45id) and (dih[k].atomids[1]==n22id)) or\
          ((dih[k].atomids[2]==n22id) and (dih[k].atomids[3]==h45id)):  # HNC48C48/C48C48NH or HNC25C13/C13C25NH or HNC25H/HC25NH 
       if dih[k].type==0:   # to prevent del a repeated dih for the same N22
        print '######## Skip deleting a Dih' 
       else:
        dih[k].type = 0
        dih1[k].type = 0
        d1x22 = d1x22 + 1                     # HNCC
        deldih22 = deldih22 + 1
        #print d1x22,'deldih HNCC  n22' 
      elif ((dih[k].atomids[0]==o1id) and (dih[k].atomids[1]==c25id)) or\
          ((dih[k].atomids[2]==c25id) and (dih[k].atomids[3]==o1id)):  # OC25CC/OC25CH or CCC25O/HCC25O
       dih[k].type = 0
       dih1[k].type = 0
       d2x22 = d2x22 + 1 
       deldih22 = deldih22 + 1                   
       #print d2x22,'deldih OCCC/OCCH n22' 
      elif ((dih[k].atomids[1]==o1id) and (dih[k].atomids[2]==c25id)) or\
          ((dih[k].atomids[1]==c25id) and (dih[k].atomids[2]==o1id)):  # HC25OC
       dih[k].type = 0
       dih1[k].type = 0
       d3x22 = d3x22 + 1  
       deldih22 = deldih22 + 1                  
       #print d3x22,'deldih HC25OC 22' 
      elif ((dih[k].atomids[0]==c25id) and (dih[k].atomids[1]==o1id) and (w[a[dih[k].atomids[3]].type]==1.008)) or\
          ((dih[k].atomids[2]==o1id) and (dih[k].atomids[3]==c25id) and (w[a[dih[k].atomids[0]].type]==1.008)):  # HCOC25/C25OCH
       dih[k].type = 0
       dih1[k].type = 0
       d4x22 = d4x22 + 1  
       deldih22 = deldih22 + 1                   
       #print d4x22,'deldih HCOC25 n22' 
      elif ((dih[k].atomids[0]==c25id) and (dih[k].atomids[1]==o1id) and (w[a[dih[k].atomids[3]].type]==12.011)) or\
          ((dih[k].atomids[2]==o1id) and (dih[k].atomids[3]==c25id) and (w[a[dih[k].atomids[0]].type]==12.011)):  # C25OCC/CCC25O
       dih[k].type = 0
       dih1[k].type = 0
       d5x22 = d5x22 + 1  
       deldih22 = deldih22 + 1                  
       #print d5x22,'deldih CCOC25  n22' 
##################################################################### XXXXXXXXXXXXX
# IMPORTANT WARNING!                                             ####
# THIS STEP IS ESSENTIAL TO DELETE REDUNDUNT(REPEATED) DIHEDRALS ####
##################################################################### 
repeateddih = 0
for i in dih1:
 for j in dih1:
  if (i != j) and (dih1[i].type==dih1[j].type) and (dih1[i].atomids[0]==dih1[j].atomids[0])\
     and (dih1[i].atomids[1]==dih1[j].atomids[1]) and (dih1[i].atomids[2]==dih1[j].atomids[2])\
     and (dih1[i].atomids[3]==dih1[j].atomids[3]):
   dih1[j].type = 0 
   repeateddih = repeateddih +1
   #print 'WARNING!!! REPEATED DIHEDRAL IS DELETED' 
  elif (i != j) and (dih1[i].type==dih1[j].type) and (dih1[i].atomids[0]==dih1[j].atomids[3])\
     and (dih1[i].atomids[1]==dih1[j].atomids[2]) and (dih1[i].atomids[2]==dih1[j].atomids[1])\
     and (dih1[i].atomids[3]==dih1[j].atomids[0]):
   dih1[j].type = 0 
   repeateddih = repeateddih +1
   #print 'WARNING!!! REPEATED DIHEDRAL IS DELETED' 
###########################################################################
# updating Dihedrals dictionaris
# rearrange/update angles dicts (dih & dih1 using dih2)
########################################################################### 
deldih = 0
kk = 0
for i in list(dih1):
 if dih1[i].type == 0:
  del dih1[i]
  deldih = deldih + 1
for j in list(dih1):
 kk = kk + 1
 dih2[kk] = dih1[j]
dih = dih2
dih1 = dih2
m.dihedrals = dih2             
m.ndihedrals = len(dih2)
###########################################################################
# step 4: 
# Re-change the type of N27 to N21
# Re-change the type of N26 to N22
###########################################################################
for i in a:
 if a[i].type==27:
  a[i].type = 21
 elif a[i].type==26:
  a[i].type = 22
###########################################################################
# step 5: output data
# Print updated Data File
###########################################################################
write.moleculefile('x_update.mol',m)
##########################################################################
####################### Checkin DATA #####################################
##########################################################################
# Check Results     ******************************************************
print 'xxx Bonds information xxx' 
print 'N21b/it=',N21b 
print 'N22b/it=',N22b 
print 'newbonds=',newbonds 
print 'delbonds=',delbonds 
print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 
print 'xxxxxxxxxxx Angles information xxxxxxx' 
print 'xxx I- when N.type 21 xxxxxxxxxxxxxxxx' 
print 'N21a=',N21a 
print 'newang21=',newang21 
print 'delang21=',delang21 
print 'nch21=',nch21 
print 'ncc21=',ncc21 
print 'C25-N21-C48=',cnc21 
print 'hnc21=',hnc21 
print 'hoc21=',hoc21 
print 'xxx II- when N.type 22 xxxxxxxxxxxxxxxx' 
print 'N22a=',N22a 
print 'newang22=',newang22 
print 'delang22=',delang22 
print 'nch22=',nch22 
print 'ncc22=',ncc22 
print 'C25-N22-C25=',cnc22a 
print 'C25-N22-C48=',cnc22b 
print 'hoc22=',hoc22 
print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 
print 'xxxxxxxxxxx Dihedrals information xxxxxxx' 
print 'xxx I- when N.type 21 xxxxxxxxxxxxxxxx' 
print 'N21d bonds=',N21d 
print 'newdih21=',newdih21 
print 'totaldeldih21=',deldih21 
print "new dih5021=",dih5021  
print "new dih5121=",dih5121 
print 'new dih5221=',dih5221 
print 'new dih5321=',dih5321 
print 'new dih5421=',dih5421 
print 'new dih5521=',dih5521 
print 'new dih5621=',dih5621 
print 'new dih5721=',dih5721 
print 'new dih5821=',dih5821 
print 'new dih5921=',dih5921 
print 'deldih21 HNCC=',d1x21 
print 'deldih21 OCCC/OCCH=',d2x21 
print 'deldih21 HC25OC=',d3x21 
print 'deldih21 HCOC25=',d4x21 
print 'deldih21 CCOC25=',d5x21 
print 'xxx I- when N.type 22 xxxxxxxxxxxxxxxx' 
print 'N22d bonds=',N22d 
print 'newdih22=',newdih22 
print 'deldih22=',deldih22 
print 'new dih5222=',dih5222 
print 'new dih5322=',dih5322 
print 'new dih5422=',dih5422 
print 'new dih5522=',dih5522 
print 'new dih5622=',dih5622 
print 'new dih5722=',dih5722 
print 'new dih5822=',dih5822 
print 'new dih5922=',dih5922 
print 'new dih6022=',dih6022 
print 'new dih6122=',dih6122 
print 'deldih22 HNCC=',d1x22 
print 'deldih22 OCCC/OCCH=',d2x22 
print 'deldih22 HC25OC=',d3x22 
print 'deldih22 HCOC25=',d4x22 
print 'deldih22 CCOC25=',d5x22 
print 'totaldeldih=',deldih 
print 'Repeatedang=',repeatedang 
print 'Repeateddih=',repeateddih 
print 'DONE!  :-) ' 
print '##########################################################################' 
print '######################### Calculating % of xlinking#######################' 
print '##########################################################################' 
k=0
bondNC24=0
bondNC25=0
for i in b:
 if (b[i].type==24):
  k=k+1
  bondNC24=bondNC24+1
 elif (b[i].type==25):
  k=k+1
  bondNC25=bondNC25+1
Xlinkpercent = k/1.28  # k is No of CN created bonds  & 1.28 is 1/128 * 100%   where 128 is the max allowed bonds
print 'No. of new CN bonds =', k 
print 'bondNC24=', bondNC24
print 'bondNC25=', bondNC25
print 'CROSSlinkpercent=', Xlinkpercent , '%'

