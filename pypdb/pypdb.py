'''
PyPDB: A Python API for the RCSB Protein Data Bank

-----

GitHub: https://github.com/williamgilpin/pypdb

PyPI: https://pypi.python.org/pypi/pypdb

-----

If you find this code useful, please consider citing the paper:

    Gilpin, W. "PyPDB: A Python API for the Protein Data Bank."
    Bioinformatics, Oxford Journals, 2015.

'''
from collections import OrderedDict, Counter
from itertools import repeat, chain
import time
import re
import json
import warnings

from pypdb.util import http_requests
from pypdb.clients.fasta import fasta_client
from pypdb.clients.pdb import pdb_client
from pypdb.clients.search import search_client
from pypdb.clients.search.operators import sequence_operators

warnings.simplefilter('always', DeprecationWarning)


# New imports needed for the updated API
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators


'''
=================
Functions for searching the RCSB PDB for lists of PDB IDs
=================
'''


class Query(object):
    """

    xThis objects takes search terms and specifications and creates object that
    can be used to query the Protein Data Bank

    Parameters
    ----------
    search_term : str

        The specific term to search in the database. For specific query types,
        the strings that will yield valid results are limited to:

        'HoldingsQuery' : A Ggeneral search of the metadata associated with PDB IDs

        'ExpTypeQuery' : Experimental Method such as 'X-RAY', 'SOLID-STATE NMR', etc

        'AdvancedKeywordQuery' : Any string that appears in the title or abstract

        'StructureIdQuery' :  Perform a search for a specific Structure ID

        'ModifiedStructuresQuery' : Search for related structures

        'AdvancedAuthorQuery' : Search by the names of authors associated with entries

        'MotifQuery' : Search for a specific motif

        'NoLigandQuery' : Find full list of PDB IDs without free ligrands

    query_type : str

        The type of query to perform, the easiest is an AdvancedKeywordQuery but more
        specific types of searches may also be performed

    return type : str
        The type of search result to return. Default "entry" returns a list of PDB IDs

    scan_params (optional) : dict()
            A dictionary containing an explicit nested search term. Use this option if you want to
            use pypdb's rate handling and other functions, but need to structure a complex JSON
            query not covered in the existing python package

    Examples
    --------

    >>> found_pdbs = Query('actin network').search()
    >>> print(found_pdbs)
    ['1D7M', '3W3D', '4A7H', '4A7L', '4A7N']

    >>> found_pdbs = Query('3W3D', query_type='ModifiedStructuresQuery').search()
    >>> print(found_pdbs[:5])
    ['1A2N', '1ACF', '1AGX', '1APM', '1ARL']

    >>> found_pdbs = found_pdbs = Query('T[AG]AGGY', query_type='MotifQuery').search()
    >>> print(found_pdbs)
    ['3LEZ', '3SGH', '4F47']

    """
    def __init__(self,
                 search_term,
                 query_type="full_text",
                 return_type="entry",
                 scan_params=None):
        """See help(Query) for documentation"""

        if query_type == "PubmedIdQuery":
            query_type = "text"
            query_subtype = "pmid"
        elif query_type == "TreeEntityQuery":
            query_type = "text"
            query_subtype = "taxid"
        elif query_type == "ExpTypeQuery":
            query_type = "text"
            query_subtype = "experiment_type"
            search_term = search_term.upper()
            if search_term not in [
                    "X-RAY DIFFRACTION", "ELECTRON MICROSCOPY",
                    "SOLID-STATE NMR", "SOLUTION NMR", "NEUTRON DIFFRACTION",
                    "ELECTRON CRYSTALLOGRAPHY", "POWDER DIFFRACTION",
                    "FIBER DIFFRACTION", "SOLUTION SCATTERING", "EPR",
                    "FLUORESCENCE TRANSFER", "INFRARED SPECTROSCOPY",
                    "THEORETICAL MODEL"
            ]:
                warnings.warn(
                    "Experimental type not recognized, search may fail .")
        elif query_type == "AdvancedAuthorQuery":
            query_type = "text"
            query_subtype = "author"
        elif query_type == "OrganismQuery":
            query_type = "text"
            query_subtype = "organism"
        elif query_type == "pfam":
            query_type = "text"
            query_subtype = "pfam"
        elif query_type == "uniprot":
            query_type = "text"
            query_subtype = "uniprot"
        else:
            query_subtype = None

        assert query_type in {
            "full_text", "text", "structure", "sequence", "seqmotif", "chemical"
        }, "Query type %s not recognized." % query_type

        assert return_type in {"entry", "polymer_entity"
                               }, "Return type %s not supported." % return_type

        self.query_type = query_type
        self.search_term = search_term
        self.return_type = return_type
        self.url = "https://search.rcsb.org/rcsbsearch/v2/query?json="
        composite_query = False
        if not scan_params:
            query_params = dict()
            query_params["type"] = "terminal"
            query_params["service"] = query_type

            if query_type in ["full_text", "text"]:
                query_params['parameters'] = {"value": search_term}

            elif query_type == "sequence":
                query_params['parameters'] = {
                    "target": "pdb_protein_sequence",
                    "value": search_term
                }
            elif query_type == "structure":
                query_params['parameters'] = {
                    "operator": "relaxed_shape_match",
                    "value": {
                        "entry_id": search_term,
                        "assembly_id": "1"
                    }
                }

#             elif query_type=='AdvancedAuthorQuery':
#                 query_params['description'] = 'Author Name: '+ search_term
#                 query_params['searchType'] = 'All Authors'
#                 query_params['audit_author.name'] = search_term
#                 query_params['exactMatch'] = 'false'

#             elif query_type=='MotifQuery':
#                 query_params['description'] = 'Motif Query For: '+ search_term
#                 query_params['motif'] = search_term

#             elif query_type=='OrganismQuery':
# #                 query_params['version'] = "B0905"
#                 query_params['description'] = 'Organism Search: Organism Name='+ search_term
#                 query_params['organismName'] = search_term
# #                 composite_query = True
#             elif query_type=='TreeEntityQuery':
#                 query_params['t'] = "1"
#                 query_params['description'] = 'TaxonomyTree Search for OTHER SEQUENCES'
#                 query_params['n'] = search_term
#                 query_params['nodeDesc'] = "OTHER SEQUENCES"

#             elif query_type in ['StructureIdQuery','ModifiedStructuresQuery']:
#                 query_params['structureIdList'] = search_term

#             elif query_type=='ExpTypeQuery':
#                 query_params['experimentalMethod'] = search_term
#                 query_params['description'] = 'Experimental Method Search : Experimental Method='+ search_term
#                 query_params['mvStructure.expMethod.value']= search_term
            if query_subtype:

                if query_subtype == "pmid":
                    query_params['parameters'] = {
                        "operator": "in",
                        "negation": False,
                        "value": [search_term],
                        "attribute":
                        "rcsb_pubmed_container_identifiers.pubmed_id"
                    }
                if query_subtype == "taxid":
                    query_params['parameters'] = {
                        "operator":
                        "exact_match",
                        "negation":
                        False,
                        "value":
                        str(search_term),
                        "attribute":
                        "rcsb_entity_source_organism.taxonomy_lineage.id"
                    }
                if query_subtype == "experiment_type":
                    query_params['parameters'] = {
                        "operator": "exact_match",
                        "negation": False,
                        "value": str(search_term),
                        "attribute": "exptl.method"
                    }
                if query_subtype == "author":
                    query_params['parameters'] = {
                        "operator": "exact_match",
                        "negation": False,
                        "value": str(search_term),
                        "attribute": "rcsb_primary_citation.rcsb_authors"
                    }
                if query_subtype == "organism":
                    query_params['parameters'] = {
                        "operator":
                        "contains_words",
                        "negation":
                        False,
                        "value":
                        str(search_term),
                        "attribute":
                        "rcsb_entity_source_organism.taxonomy_lineage.name"
                    }
                if query_subtype == "pfam":
                    query_params['parameters'] = {
                        "operator": "exact_match",
                        "negation": False,
                        "value": str(search_term),
                        "attribute":
                        "rcsb_polymer_entity_annotation.annotation_id"
                    }
                if query_subtype == "uniprot":
                    query_params['parameters'] = {
                        "operator": "exact_match",
                        "negation": False,
                        "value": str(search_term),
                        "attribute": 
                        "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
                        }

            self.scan_params = dict()
            self.scan_params["query"] = query_params
            self.scan_params["return_type"] = return_type
            self.scan_params["request_options"] = {"results_verbosity": "verbose"} # v2

            if return_type == "entry":
                self.scan_params["request_options"] = {"return_all_hits": True}

        else:
            self.scan_params = scan_params

    def search(self, num_attempts=1, sleep_time=0.5):
        """
        Perform a search of the Protein Data Bank using the REST API

        Parameters
        ----------

        num_attempts : int
            In case of a failed retrieval, the number of attempts to try again
        sleep_time : int
            The amount of time to wait between requests, in case of
            API rate limits
        """

        query_text = json.dumps(self.scan_params, indent=4)
        response = http_requests.request_limited(self.url,
                                                 rtype="POST",
                                                 data=query_text)

        if response is None or response.status_code != 200:
            warnings.warn("Retrieval failed, returning None")
            return None

        response_val = json.loads(response.text)

        if self.return_type == "entry":
            idlist = walk_nested_dict(response_val,
                                      "identifier",
                                      maxdepth=25,
                                      outputs=[])
            return idlist
        else:
            return response_val


# def do_search(scan_params):
#     '''Convert dict() to XML object an then send query to the RCSB PDB

#     This function takes a valid query dict() object, converts it to XML,
#     and then sends a request to the PDB for a list of IDs corresponding
#     to search results

#     Parameters
#     ----------

#     scan_params : dict
#         A dictionary of query attributes to use for
#         the search of the PDB

#     Returns
#     -------

#     idlist : list
#         A list of PDB ids returned by the search

#     Examples
#     --------
#     This method usually gets used in tandem with make_query

#     >>> a = make_query('actin network')
#     >>> print (a)
#     {'orgPdbQuery': {'description': 'Text Search for: actin',
#     'keywords': 'actin',
#     'queryType': 'AdvancedKeywordQuery'}}

#     >>> search_dict = make_query('actin network')
#     >>> found_pdbs = do_search(search_dict)
#     >>> print(found_pdbs)
#     ['1D7M', '3W3D', '4A7H', '4A7L', '4A7N']

#     >>> search_dict = make_query('T[AG]AGGY',querytype='MotifQuery')
#     >>> found_pdbs = do_search(search_dict)
#     >>> print(found_pdbs)
#     ['3LEZ', '3SGH', '4F47']
#     '''
#     q = Query('search_term', 'HoldingsQuery', scan_params=scan_params)
#     return q.search()

# def do_protsym_search(point_group, min_rmsd=0.0, max_rmsd=7.0):
#     '''Performs a protein symmetry search of the PDB

#     This function can search the Protein Data Bank based on how closely entries
#     match the user-specified symmetry group

#     Parameters
#     ----------

#     point_group : str
#         The name of the symmetry point group to search. This includes all the standard
#         abbreviations for symmetry point groups (e.g., C1, C2, D2, T, O, I, H, A1)

#     min_rmsd : float
#         The smallest allowed total deviation (in Angstroms) for a result
#         to be classified as having a matching symmetry

#     max_rmsd : float
#         The largest allowed total deviation (in Angstroms) for a result
#         to be classified as having a matching symmetry

#     Returns
#     -------

#     idlist : list of strings
#         A list of PDB IDs resulting from the search

#     Examples
#     --------

#     >>> kk = do_protsym_search('C9', min_rmsd=0.0, max_rmsd=1.0)
#     >>> print(kk[:5])
#     ['1KZU', '1NKZ', '2FKW', '3B8M', '3B8N']

#     '''
#     query_params = dict()
#     query_params['queryType'] = 'PointGroupQuery'
#     query_params['rMSDComparator'] = 'between'

#     query_params['pointGroup'] = point_group
#     query_params['rMSDMin'] = min_rmsd
#     query_params['rMSDMax'] = max_rmsd

#     scan_params = dict()
#     scan_params['orgPdbQuery'] = query_params
#     idlist =  do_search(scan_params)
#     return idlist

# def get_all():
#     """Return a list of all PDB entries currently in the RCSB Protein Data Bank

#     Returns
#     -------

#     out : list of str
#         A list of all of the PDB IDs currently in the RCSB PDB

#     Examples
#     --------

#     >>> print(get_all()[:10])
#     ['100D', '101D', '101M', '102D', '102L', '102M', '103D', '103L', '103M', '104D']

#     """

#     url = 'http://www.rcsb.org/pdb/rest/getCurrent'
#     #response = requests.get(url)
#     response = http_requests.request_limited(url)

#     if response.status_code == 200:
#         pass
#     else:
#         warnings.warn("Retrieval failed, returning None")
#         return None

#     result  = str(response.text)

#     p = re.compile('structureId=\"...."')
#     matches = p.findall(str(result))
#     out = list()
#     for item in matches:
#         out.append(item[-5:-1])

#     return out
'''
=================
Functions for looking up information given PDB ID
=================
'''


def get_info(pdb_id, url_root='https://data.rcsb.org/rest/v1/core/entry/'):
    '''Look up all information about a given PDB ID

    Parameters
    ----------

    pdb_id : string
        A 4 character string giving a pdb entry of interest

    url_root : string
        The string root of the specific url for the request type

    Returns
    -------

    out : dict()
        An ordered dictionary object corresponding to entry information

    '''
    pdb_id = pdb_id.replace(":", "/")  # replace old entry identifier
    url = url_root + pdb_id
    response = http_requests.request_limited(url)

    if response is None or response.status_code != 200:
        warnings.warn("Retrieval failed, returning None")
        return None

    result = str(response.text)

    out = json.loads(result)

    return out


get_all_info = get_info  # Alias
describe_pdb = get_info  # Alias for now; eventually make this point to the Graph search https://data.rcsb.org/migration-guide.html#pdb-file-description
get_entity_info = get_info  # Alias


def get_pdb_file(pdb_id: str, filetype='pdb', compression=False):
    """Deprecated wrapper for fetching PDB files from RCSB Database.

    For new uses, please use `pypdb/clients/pdb/pdb_client.py`
    """

    warnings.warn(
        "The `get_pdb_file` function within pypdb.py is deprecated."
        "See `pypdb/clients/pdb/pdb_client.py` for a near-identical "
        "function to use", DeprecationWarning)

    if filetype == 'pdb':
        filetype_enum = pdb_client.PDBFileType.PDB
    elif filetype == 'cif':
        filetype_enum = pdb_client.PDBFileType.CIF
    elif filetype == 'xml':
        filetype_enum = pdb_client.PDBFileType.XML
    elif filetype == 'structfact':
        filetype_enum = pdb_client.PDBFileType.STRUCTFACT
    else:
        warnings.warn(
            "Filetype specified to `get_pdb_file` appears to be invalid")

    return pdb_client.get_pdb_file(pdb_id, filetype_enum, compression)


# https://data.rcsb.org/migration-guide.html#chem-comp-description
def describe_chemical(chem_id):
#     """

#     Parameters
#     ----------

#     chem_id : string
#         A 3 character string representing the full chemical sequence of interest (ie, NAG)

#     Returns
#     -------

#     out : dict
#         A dictionary containing the chemical description associated with the PDB ID

#     Examples
#     --------
#     >>> chem_desc = describe_chemical('NAG')
#     >>> print(chem_desc["rcsb_chem_comp_descriptor"]["smiles"])
#     'CC(=O)NC1C(C(C(OC1O)CO)O)O'
#     """
    if (len(chem_id) > 3):
        raise Exception("Ligand id with more than 3 characters provided")

    return get_info(chem_id, url_root = 'https://data.rcsb.org/rest/v1/core/chemcomp/')

# def get_ligands(pdb_id):
#     """Return ligands of given PDB ID

#     Parameters
#     ----------

#     pdb_id : string
#         A 4 character string giving a pdb entry of interest

#     Returns
#     -------

#     out : dict
#         A dictionary containing a list of ligands associated with the entry

#     Examples
#     --------
#     >>> ligand_dict = get_ligands('100D')
#     >>> print(ligand_dict)
#     {'id': '100D',
#     'ligandInfo': {'ligand': {'@chemicalID': 'SPM',
#                            '@molecularWeight': '202.34',
#                            '@structureId': '100D',
#                            '@type': 'non-polymer',
#                            'InChI': 'InChI=1S/C10H26N4/c11-5-3-9-13-7-1-2-8-14-10-4-6-12/h13-14H,1-12H2',
#                            'InChIKey': 'PFNFFQXMRSDOHW-UHFFFAOYSA-N',
#                            'chemicalName': 'SPERMINE',
#                            'formula': 'C10 H26 N4',
#                            'smiles': 'C(CCNCCCN)CNCCCN'}}}

#     """
#     out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/ligandInfo?structureId=')
#     out = to_dict(out)
#     return remove_at_sign(out['structureId'])

# def get_gene_onto(pdb_id):
#     """Return ligands of given PDB_ID

#     Parameters
#     ----------

#     pdb_id : string
#         A 4 character string giving a pdb entry of interest

#     Returns
#     -------

#     out : dict
#         A dictionary containing the gene ontology information associated with the entry

#     Examples
#     --------

#     >>> gene_info = get_gene_onto('4Z0L')
#     >>> print(gene_info['term'][0])
#     {'@chainId': 'A',
#      '@id': 'GO:0001516',
#      '@structureId': '4Z0L',
#      'detail': {'@definition': 'The chemical reactions and pathways resulting '
#                                'in the formation of prostaglandins, any of a '
#                                'group of biologically active metabolites which '
#                                'contain a cyclopentane ring.',
#                 '@name': 'prostaglandin biosynthetic process',
#                 '@ontology': 'B',
#                 '@synonyms': 'prostaglandin anabolism, prostaglandin '
#                              'biosynthesis, prostaglandin formation, '
#                              'prostaglandin synthesis'}}
#     """
#     out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/goTerms?structureId=')
#     out = to_dict(out)
#     if not out['goTerms']:
#         return None
#     out = remove_at_sign(out['goTerms'])
#     return out

# def get_seq_cluster(pdb_id_chain):
#     """Get the sequence cluster of a PDB ID plus a pdb_id plus a chain,

#     Parameters
#     ----------

#     pdb_id_chain : string
#         A string denoting a 4 character PDB ID plus a one character chain
#         offset with a dot: XXXX.X, as in 2F5N.A

#     Returns
#     -------

#     out : dict
#         A dictionary containing the sequence cluster associated with the PDB
#         entry and chain

#     Examples
#     --------

#     >>> sclust = get_seq_cluster('2F5N.A')
#     >>> print(sclust['pdbChain'][:10])
#     [{'@name': '4PD2.A', '@rank': '1'},
#      {'@name': '3U6P.A', '@rank': '2'},
#      {'@name': '4PCZ.A', '@rank': '3'},
#      {'@name': '3GPU.A', '@rank': '4'},
#      {'@name': '3JR5.A', '@rank': '5'},
#      {'@name': '3SAU.A', '@rank': '6'},
#      {'@name': '3GQ4.A', '@rank': '7'},
#      {'@name': '1R2Z.A', '@rank': '8'},
#      {'@name': '3U6E.A', '@rank': '9'},
#      {'@name': '2XZF.A', '@rank': '10'}]

#     """

#     url_root = 'http://www.rcsb.org/pdb/rest/sequenceCluster?structureId='
#     out = get_info(pdb_id_chain, url_root = url_root)
#     out = to_dict(out)
#     return remove_at_sign(out['sequenceCluster'])


def get_blast(pdb_id, chain_id='A', identity_cutoff=0.99, verbosity=True):
    """
    ---
    WARNING: this function is deprecated and slated to be deleted due to RCSB
    API changes.

    See `pypdb/clients/search/EXAMPLES.md` for examples to use a
    `SequenceOperator` search to similar effect
    ---

    Return BLAST search results for a given PDB ID.

    Parameters
    ----------
    pdb_id : string
        A 4 character string giving a pdb entry of interest

    chain_id : string
        A single character designating the chain ID of interest
    identity_cutoff: float
        Identity % at which to cut off results.


    Returns
    -------

    out : List of PDB IDs that match the given search.

    Examples
    --------

    >>> blast_results = get_blast('2F5N', chain_id='A')
    >>> print(blast_results[50])
    PELPEVETVRRELEKRIVGQKIISIEATYPRMVL--TGFEQLKKELTGKTIQGISRRGKYLIFEIGDDFRLISHLRMEGKYRLATLDAPREKHDHL
    TMKFADG-QLIYADVRKFGTWELISTDQVLPYFLKKKIGPEPTYEDFDEKLFREKLRKSTKKIKPYLLEQTLVAGLGNIYVDEVLWLAKIHPEKET
    NQLIESSIHLLHDSIIEILQKAIKLGGSSIRTY-SALGSTGKMQNELQVYGKTGEKCSRCGAEIQKIKVAGRGTHFCPVCQQ
    """

    warnings.warn(
        "The `get_blast` function is slated for deprecation."
        "See `pypdb/clients/search/EXAMPLES.md` for examples to use a"
        "`SequenceOperator` search to similar effect", DeprecationWarning)

    fasta_entries = fasta_client.get_fasta_from_rcsb_entry(pdb_id)
    valid_sequences = [
        fasta_entry.sequence for fasta_entry in fasta_entries
        if chain_id in fasta_entry.chains
    ]

    matches_any_sequence_in_chain_query = search_client.QueryGroup(
        logical_operator=search_client.LogicalOperator.OR, queries=[])
    for valid_sequence in valid_sequences:
        matches_any_sequence_in_chain_query.queries.append(
            sequence_operators.SequenceOperator(
                sequence=valid_sequence,
                identity_cutoff=identity_cutoff,
                evalue_cutoff=1000))

    return search_client.perform_search_with_graph(
        query_object=matches_any_sequence_in_chain_query,
        return_raw_json_dict=True)


# def get_pfam(pdb_id):
#     """Return PFAM annotations of given PDB_ID

#     Parameters
#     ----------

#     pdb_id : string
#         A 4 character string giving a pdb entry of interest

#     Returns
#     -------

#     out : dict
#         A dictionary containing the PFAM annotations for the specified PDB ID

#     Examples
#     --------

#     >>> pfam_info = get_pfam('2LME')
#     >>> print(pfam_info)
#     {'pfamHit': {'@pfamAcc': 'PF03895.10', '@pfamName': 'YadA_anchor',
#     '@structureId': '2LME', '@pdbResNumEnd': '105', '@pdbResNumStart': '28',
#     '@pfamDesc': 'YadA-like C-terminal region', '@eValue': '5.0E-22', '@chainId': 'A'}}

#     """
#     out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/hmmer?structureId=')
#     out = to_dict(out)
#     if not out['hmmer3']:
#         return dict()
#     return remove_at_sign(out['hmmer3'])

# def get_clusters(pdb_id):
#     """Return cluster related web services of given PDB_ID

#     Parameters
#     ----------

#     pdb_id : string
#         A 4 character string giving a pdb entry of interest

#     Returns
#     -------

#     out : dict
#         A dictionary containing the representative clusters for the specified PDB ID

#     Examples
#     --------

#     >>> clusts = get_clusters('4hhb.A')
#     >>> print(clusts)
#     {'pdbChain': {'@name': '2W72.A'}}

#     """
#     out = get_info(pdb_id, url_root = 'http://www.rcsb.org/pdb/rest/representatives?structureId=')
#     out = to_dict(out)
#     return remove_at_sign(out['representatives'])


def find_results_gen(search_term, field='title'):
    '''
    Return a generator of the results returned by a search of
    the protein data bank. This generator is used internally.

    Parameters
    ----------

    search_term : str
        The search keyword

    field : str
        The type of information to record about each entry

    Examples
    --------

    >>> result_gen = find_results_gen('bleb')
    >>> pprint.pprint([item for item in result_gen][:5])
    ['MYOSIN II DICTYOSTELIUM DISCOIDEUM MOTOR DOMAIN S456Y BOUND WITH MGADP-BEFX',
     'MYOSIN II DICTYOSTELIUM DISCOIDEUM MOTOR DOMAIN S456Y BOUND WITH MGADP-ALF4',
     'DICTYOSTELIUM DISCOIDEUM MYOSIN II MOTOR DOMAIN S456E WITH BOUND MGADP-BEFX',
     'MYOSIN II DICTYOSTELIUM DISCOIDEUM MOTOR DOMAIN S456E BOUND WITH MGADP-ALF4',
     'The structural basis of blebbistatin inhibition and specificity for myosin '
     'II']

    '''
    search_result_ids = Query(search_term).search()

    all_titles = []
    for pdb_id in search_result_ids:
        result = get_info(pdb_id)
        if field in result.keys():
            yield result[field]


def find_papers(search_term, max_results=10, **kwargs):
    '''
    Return an ordered list of the top papers returned by a keyword search of
    the RCSB PDB

    Parameters
    ----------

    search_term : str
        The search keyword

    max_results : int
        The maximum number of results to return

    Returns
    -------

    all_papers : list of strings
        A descending-order list containing the top papers associated with
        the search term in the PDB

    Examples
    --------

    >>> matching_papers = find_papers('crispr',max_results=3)
    >>> print(matching_papers)
    ['Crystal structure of a CRISPR-associated protein from thermus thermophilus',
    'CRYSTAL STRUCTURE OF HYPOTHETICAL PROTEIN SSO1404 FROM SULFOLOBUS SOLFATARICUS P2',
    'NMR solution structure of a CRISPR repeat binding protein']

    '''
    all_papers = list()
    id_list = Query(search_term).search()
    for pdb_id in id_list[:max_results]:
        pdb_info = get_info(pdb_id)
        all_papers += [item["title"] for item in pdb_info["citation"]]
    return remove_dupes(all_papers)


# def find_authors(search_term, **kwargs):
#     '''Return an ordered list of the top authors returned by a keyword search of
#     the RCSB PDB

#     This function is based on the number of unique PDB entries a given author has
#     his or her name associated with, and not author order or the ranking of the
#     entry in the keyword search results. So if an author tends to publish on topics
#     related to the search_term a lot, even if those papers are not the best match for
#     the exact search, he or she will have priority in this function over an author
#     who wrote the one paper that is most relevant to the search term. For the latter
#     option, just do a standard keyword search using do_search.

#     Parameters
#     ----------

#     search_term : str
#         The search keyword

#     max_results : int
#         The maximum number of results to return

#     Returns
#     -------

#     out : list of str

#     Examples
#     --------

#     >>> top_authors = find_authors('crispr',max_results=100)
#     >>> print(top_authors[:10])
#     ['Doudna, J.A.', 'Jinek, M.', 'Ke, A.', 'Li, H.', 'Nam, K.H.']

#     '''

#     all_individuals = parse_results_gen(search_term, field='citation_authors', **kwargs)

#     full_author_list = []
#     for individual in all_individuals:
#         individual = individual.replace('.,', '.;')
#         author_list_clean = [x.strip() for x in individual.split(';')]
#         full_author_list+=author_list_clean

#     out = list(chain.from_iterable(repeat(ii, c) for ii,c in Counter(full_author_list).most_common()))

#     return remove_dupes(out)

# def find_dates(search_term, **kwargs):
#     '''
#     Return an ordered list of the PDB submission dates returned by a
#     keyword search of the RCSB PDB. This can be used to assess the
#     popularity of a gievne keyword or topic

#     Parameters
#     ----------

#     search_term : str
#         The search keyword

#     max_results : int
#         The maximum number of results to return

#     Returns
#     -------

#     all_dates : list of str
#         A list of calendar strings associated with the search term, these can
#         be converted directly into time or datetime objects

#     '''
#     all_dates = parse_results_gen(search_term, field='deposition_date', **kwargs)
#     return all_dates
'''
=================
Helper Functions
=================
'''


def to_dict(odict):
    '''Convert OrderedDict to dict

    Takes a nested, OrderedDict() object and outputs a
    normal dictionary of the lowest-level key:val pairs

    Parameters
    ----------

    odict : OrderedDict

    Returns
    -------

    out : dict

        A dictionary corresponding to the flattened form of
        the input OrderedDict

    '''

    out = json.loads(json.dumps(odict))
    return out


def remove_at_sign(kk):
    '''Remove the '@' character from the beginning of key names in a dict()

    Parameters
    ----------

    kk : dict
        A dictionary containing keys with the @ character
        (this pops up a lot in converted XML)

    Returns
    -------

    kk : dict (modified in place)
        A dictionary where the @ character has been removed

    '''
    tagged_keys = [thing for thing in kk.keys() if thing.startswith('@')]
    for tag_key in tagged_keys:
        kk[tag_key[1:]] = kk.pop(tag_key)

    return kk


def remove_dupes(list_with_dupes):
    '''Remove duplicate entries from a list while preserving order

    This function uses Python's standard equivalence testing methods in
    order to determine if two elements of a list are identical. So if in the list [a,b,c]
    the condition a == b is True, then regardless of whether a and b are strings, ints,
    or other, then b will be removed from the list: [a, c]

    Parameters
    ----------

    list_with_dupes : list
        A list containing duplicate elements

    Returns
    -------
    out : list
        The list with the duplicate entries removed by the order preserved


    Examples
    --------
    >>> a = [1,3,2,4,2]
    >>> print(remove_dupes(a))
    [1,3,2,4]

    '''
    visited = set()
    visited_add = visited.add
    out = [
        entry for entry in list_with_dupes
        if not (entry in visited or visited_add(entry))
    ]
    return out


def walk_nested_dict(my_result, term, outputs=[], depth=0, maxdepth=25):
    '''
    For a nested dictionary that may itself comprise lists of
    dictionaries of unknown length, determine if a key is anywhere
    in any of the dictionaries using a depth-first search

    Parameters
    ----------

    my_result : dict
        A nested dict containing lists, dicts, and other objects as vals

    term : str
        The name of the key stored somewhere in the tree

    maxdepth : int
        The maximum depth to search the results tree

    depth : int
        The depth of the search so far.
        Users don't usually access this.

    outputs : list
        All of the positive search results collected so far.
        Users don't usually access this.

    Returns
    -------

    outputs : list
        All of the search results.

    '''

    if depth > maxdepth:
        warnings.warn(
            'Maximum recursion depth exceeded. Returned None for the search results,'
            + ' try increasing the maxdepth keyword argument.')
        return None

    depth = depth + 1

    if type(my_result) == dict:
        if term in my_result.keys():
            outputs.append(my_result[term])

        else:
            new_results = list(my_result.values())
            walk_nested_dict(new_results,
                             term,
                             outputs=outputs,
                             depth=depth,
                             maxdepth=maxdepth)

    elif type(my_result) == list:
        for item in my_result:
            walk_nested_dict(item,
                             term,
                             outputs=outputs,
                             depth=depth,
                             maxdepth=maxdepth)

    else:
        pass
        # dead leaf

    # this conditional may not be necessary
    if outputs:
        return outputs
    else:
        return None
