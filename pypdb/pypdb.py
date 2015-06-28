'''
pypdb
-----
William Gilpin, 2015
Please heed the license accompanying this file.
'''

from matplotlib.pyplot import *
from numpy import *

import urllib.request
import dicttoxml
import xmltodict
import re

from json import loads, dumps

from collections import OrderedDict

# functions for conducting searches and obtaining lists of PDB ids
def make_query(search_term, querytype='AdvancedKeywordQuery'):
    '''
    This function takes a list of search terms and specifications
    and repackages it as a dictionary object to be used to conduct a search

    Inputs
    ------
    search_term : str 
    
        The specific term to search in the database. For specific query types,
        the toekns that will yeild valid results is limited. For example, a ExpTypeQuery
        will only work if the search_term is 'X-RAY', 'SOLID-STATE NMR', etc
    
    querytype : str
    
        The type of query to perform, the easiest is an AdvancedKeywordQuery but more
        specific types of searches may also be performed
        
    Returns
    -------
    
    scan_params : dict
        
        A dictionary representing the query
        
    William Gilpin, 2015. Please heed attribution described in the license
    '''
    assert querytype in {'HoldingsQuery', 'ExpTypeQuery',
                         'AdvancedKeywordQuery','StructureIdQuery',
                         'ModifiedStructuresQuery'
                        }
  
    query_params = dict()
    query_params['queryType'] = querytype
    
    if querytype=='AdvancedKeywordQuery':
        query_params['description'] = 'Text Search for: '+ search_term
        query_params['keywords'] = search_term
        
    # search for a specific structure
    elif querytype in ['StructureIdQuery','ModifiedStructuresQuery']:
        query_params['structureIdList'] = search_term
    
        
    elif querytype=='ExpTypeQuery':
        query_params['experimentalMethod'] = search_term
        query_params['description'] = 'Experimental Method Search : Experimental Method='+ search_term
        query_params['mvStructure.expMethod.value']= search_term
        
    
    scan_params = dict()
    scan_params['orgPdbQuery'] = query_params
    
    
    return scan_params

def do_search(scan_params):
    '''
    This function takes a valid query dict() object, converts it to XML,
    and then sends a request to the PDB for a list of IDs corresponding to search results
    
    Inputs
    ------
    
    scan_params : dict
        A dictionary of query attributes to use for
        the search of the PDB
    
    
    Returns
    -------
    
    idlist : list
        A list of PDB ids returned by the search
        
    William Gilpin, 2015. Please heed attribution described in the license
    '''
    
    url = 'http://www.rcsb.org/pdb/rest/search'


    queryText = dicttoxml.dicttoxml(scan_params)

    req = urllib.request.Request(url, data=queryText)
    f = urllib.request.urlopen(req)
    result = f.read()
    assert result

    idlist = str(result)
    idlist =idlist.split('\\n')
    idlist[0] = idlist[0][-4:]
    kk = idlist.pop(-1)
    
    return idlist




# Functions for cleaning stuff up
def to_dict(odict):
    '''
    Takes a nested, OrderedDict() object and returns a 
    normal dictionary of the lowest-level key:val pairs
    '''
    return loads(dumps(odict))

def remove_tags(kk):
    '''
    Removes the '@' character from the beginning of key names in a dict()
    '''
    tagged_keys = [thing for thing in kk.keys() if thing.startswith('@')]
    for tag_key in tagged_keys:
        kk[tag_key[1:]] = kk.pop(tag_key)
        
    return kk

# functions for obtaining information about PDB id files
def get_all():
    """
    Return a list of all PDB entried currently in the database
    
    """
    
    url = 'http://www.rcsb.org/pdb/rest/getCurrent'
    
    req = urllib.request.Request(url)
    f = urllib.request.urlopen(req)
    result = f.read()
    assert result
    
    kk = str(result)

    p = re.compile('structureId=\"...."')
    matches = p.findall(str(result))
    out = list()
    for item in matches:
        out.append(item[-5:-1])

    return out

def get_info(pdb_id, url_root='http://www.rcsb.org/pdb/rest/describeMol?structureId='):
    '''
    Look up all information about a given PDB ID
    
    pdb_id : string
        A 4 character string giving a pdb entry of interest
        
    url_root : string
        The string root of the specific url for the request type
        
    Returns
    -------
    
    out : OrderedDict
        An ordered dictionary object corresponding to bare xml
        
    '''
    
    url = url_root + pdb_id
    req = urllib.request.Request(url)
    f = urllib.request.urlopen(req)
    result = f.read()
    assert result
    
    out = xmltodict.parse(result,process_namespaces=True)
    
    return out
    
def describe_pdb(pdb_id):
    """
    Return a text pdb description from PDB
    """
    out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/describePDB?structureId=')
    out = to_dict(out)
    return remove_tags(out['PDBdescription']['PDB'])

def get_entity_info(pdb_id):
    """
    Return pdb id information
    """
    out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/getEntityInfo?structureId=')
    out = to_dict(out)
    return remove_tags( out['entityInfo']['PDB'] )

def describe_chemical(pdb_id):
    """
    Return pdb id information
    """
    out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/describeHet?chemicalID=')
    out = to_dict(out)
    return out

def get_blast(pdb_id):
    """
    Return BLAST of given PDB_ID
    """
    out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/getBlastPDB1?structureId=')
    out = to_dict(out)
    return out

def get_ligands(pdb_id):
    """
    Return ligands of given PDB_ID
    """
    out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/ligandInfo?structureId=')
    out = to_dict(out)
    return remove_tags(out['structureId'])

def get_gene_onto(pdb_id):
    """
    Return ligands of given PDB_ID
    """
    out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/goTerms?structureId=')
    out = to_dict(out)
    return remove_tags(out['goTerms'])

def get_pfam(pdb_id):
    """
    Return PFAM annotations of given PDB_ID
    """
    out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/hmmer?structureId=')
    out = to_dict(out)
    return remove_tags(out['hmmer3'])

def get_clusters(pdb_id):
    """
    Return cluster related web services of given PDB_ID
    """
    out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/representatives?structureId=')
    out = to_dict(out)
    return remove_tags(out['representatives'])



