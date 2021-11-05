#Geobacillus spec interface to BuildOrg2, currenly based on
#ModelFromXL, likely to change in future.
#basic structure of this program similar to MakeBarleyModel.py

from Util import Set
from Bioinf import PyoCyc

import BuildOrg2, CompartmentDic
reloads = [BuildOrg2,CompartmentDic]
#enq = lambda x: x.join(('"','"'))

org = 'Reutro' #note version db - 21.0
glo = 'MetaCyc_18.5' #note version db



comp_dict = CompartmentDic.CompartmentDic()

for r in reloads:
    reload(r)

    
def Init(o_db = None):
    global orgdb
    if o_db == None:
       orgdb = PyoCyc.Organism(data='data', path = '/home/nicole/Desktop/db/reut381666cyc/21.0/')
    else:
        orgdb = o_db
    
def GetReacs(o_db):
    rv = o_db.Reaction.keys()
    return Set.MakeSet(rv) # MakeSet to eliminate duplicates

def BuildModel(Pathways=[], db_o = None, OtherReacs=[],fname='AutoReutro21_2.spy'):
    Init(db_o)
    reacs = GetReacs(orgdb)
    BuildOrg2.BuildModel(reacs, comp_dict, orgdb, fname)

def RebuildModel(m,Pathways=[], db_o = None, OtherReacs=[], fname='AutoReutro21_2.spy'):
    Init(db_o)
    reacs = GetReacs(orgdb)
    BuildOrg2.RebuildModel(m,reacs, comp_dict, orgdb, fname)
