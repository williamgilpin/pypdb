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
import json
import warnings

import requests

from pypdb.util import http_requests
from pypdb.clients.fasta import fasta_client
from pypdb.clients.pdb import pdb_client
from pypdb.clients.search import search_client
from pypdb.clients.search.operators import chemical_operators
from pypdb.clients.search.operators import selection_operators
from pypdb.clients.search.operators import seqmotif_operators
from pypdb.clients.search.operators import sequence_operators
from pypdb.clients.search.operators import structure_operators

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

    This object takes search terms and specifications and creates an object
    that can be used to query the Protein Data Bank.

    Parameters
    ----------
    search_term : str

        The specific term to search in the database. Its interpretation
        depends on `query_type`.

    query_type : str

        The type of query to perform. The easiest is the default `full_text`
        search, but more specific types of searches may also be performed:

        'full_text' : Any string appearing anywhere in the entry's metadata

        'text' : A search against a specific structured attribute

        'sequence' : Find entries with a similar protein sequence (MMseqs2)

        'seqmotif' : Search for a specific sequence motif

        'structure' : Find entries with a similar 3D structure

        'chemical' : Search by ligand. Descriptors longer than three
            characters are treated as SMILES/InChI strings; shorter terms are
            treated as chemical component IDs (e.g. 'NAG')

        'PubmedIdQuery' : Search by the PubMed ID of the associated paper

        'TreeEntityQuery' : Search by source organism NCBI TaxId

        'ExpTypeQuery' : Experimental method, e.g. 'X-RAY DIFFRACTION'

        'AdvancedAuthorQuery' : Search by the names of associated authors

        'OrganismQuery' : Search by the name of the source organism

        'pfam' : Search by Pfam accession number, e.g. 'PF00008'

        'uniprot' : Search by UniProt accession, e.g. 'A0A023GPI8'

        'ec_number' : Search by Enzyme Classification number, e.g. '2.3.2.5'.
            Partial EC numbers match everything below them in the hierarchy,
            so '2.3.2' returns all of its sub-classes

    return_type : str
        The type of search result to return. Default "entry" returns a list of
        PDB IDs; "polymer_entity" returns the raw JSON response.

    scan_params (optional) : dict()
            A dictionary containing an explicit nested search term. Use this option if you want to
            use pypdb's rate handling and other functions, but need to structure a complex JSON
            query not covered in the existing python package

    Examples
    --------

    >>> found_pdbs = Query('actin network').search()
    >>> print(found_pdbs[:5])
    ['1D7M', '3W3D', '4A7H', '4A7L', '4A7N']

    >>> found_pdbs = Query('6239', query_type='TreeEntityQuery').search()
    >>> print(found_pdbs[:5])
    ['1D4X', '1EIL', '1EL0', '1ELE', '1EMA']

    """
    # Maps the legacy `query_type` aliases onto the "text" search service,
    # along with the query subtype used to pick an attribute below.
    _QUERY_TYPE_ALIASES = {
        "PubmedIdQuery": "pmid",
        "TreeEntityQuery": "taxid",
        "ExpTypeQuery": "experiment_type",
        "AdvancedAuthorQuery": "author",
        "OrganismQuery": "organism",
        "pfam": "pfam",
        "uniprot": "uniprot",
        "ec_number": "ec_number",
    }

    # Attribute (and matching text operator) backing each query subtype.
    _SUBTYPE_ATTRIBUTES = {
        "pmid": ("in", "rcsb_pubmed_container_identifiers.pubmed_id"),
        "taxid":
        ("exact_match", "rcsb_entity_source_organism.taxonomy_lineage.id"),
        "experiment_type": ("exact_match", "exptl.method"),
        "author": ("exact_match", "rcsb_primary_citation.rcsb_authors"),
        "organism": ("contains_words",
                     "rcsb_entity_source_organism.taxonomy_lineage.name"),
        "pfam":
        ("exact_match", "rcsb_polymer_entity_annotation.annotation_id"),
        "uniprot": (
            "exact_match",
            "rcsb_polymer_entity_container_identifiers."
            "reference_sequence_identifiers.database_accession",
        ),
        # Chemical component (ligand) IDs, e.g. "NAG"
        "chem_comp_id":
        ("exact_match",
         "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id"),
        # Enzyme Classification numbers, e.g. "2.3.2.5". This is a lineage
        # attribute, so a partial EC number such as "2.3.2" matches every
        # entry beneath it in the EC hierarchy.
        "ec_number": ("in", "rcsb_polymer_entity.rcsb_ec_lineage.id"),
    }

    EXPERIMENTAL_METHODS = [
        "X-RAY DIFFRACTION", "ELECTRON MICROSCOPY", "SOLID-STATE NMR",
        "SOLUTION NMR", "NEUTRON DIFFRACTION", "ELECTRON CRYSTALLOGRAPHY",
        "POWDER DIFFRACTION", "FIBER DIFFRACTION", "SOLUTION SCATTERING",
        "EPR", "FLUORESCENCE TRANSFER", "INFRARED SPECTROSCOPY",
        "THEORETICAL MODEL"
    ]

    def __init__(self,
                 search_term,
                 query_type="full_text",
                 return_type="entry",
                 scan_params=None):
        """See help(Query) for documentation"""
        self._uses_custom_scan_params = scan_params is not None

        query_subtype = self._QUERY_TYPE_ALIASES.get(query_type)
        if query_subtype is not None:
            if query_subtype == "experiment_type":
                search_term = search_term.upper()
                if search_term not in self.EXPERIMENTAL_METHODS:
                    warnings.warn(
                        "Experimental type not recognized, search may fail .")
            query_type = "text"

        assert query_type in {
            "full_text", "text", "structure", "sequence", "seqmotif", "chemical"
        }, "Query type %s not recognized." % query_type

        assert return_type in {"entry", "polymer_entity"
                               }, "Return type %s not supported." % return_type

        # A short chemical search term is a chemical component (ligand) ID
        # rather than a SMILES/InChI descriptor, and is served by a text search
        # over the component ID attribute.
        if query_type == "chemical" and not self._is_chemical_descriptor(
                search_term):
            query_type = "text"
            query_subtype = "chem_comp_id"
            search_term = str(search_term).upper()

        self.query_type = query_type
        self.search_term = search_term
        self.return_type = return_type
        self.query_subtype = query_subtype
        self.url = "https://search.rcsb.org/rcsbsearch/v2/query?json="
        self.scan_params = scan_params if scan_params else None

    @staticmethod
    def _is_chemical_descriptor(search_term):
        """Whether a chemical search term is a SMILES/InChI descriptor.

        Chemical component IDs in the PDB are at most three characters long,
        so longer terms are treated as structural descriptors.
        """
        search_term = str(search_term)
        return search_term.startswith("InChI=") or len(search_term) > 3

    def _get_return_type_enum(self) -> search_client.ReturnType:
        if self.return_type == "entry":
            return search_client.ReturnType.ENTRY
        if self.return_type == "polymer_entity":
            return search_client.ReturnType.POLYMER_ENTITY
        raise AssertionError("Return type %s not supported." % self.return_type)

    def _build_search_operator(self):
        if self.query_type == "full_text":
            return text_operators.DefaultOperator(value=str(self.search_term))

        if self.query_type == "sequence":
            return sequence_operators.SequenceOperator(
                sequence=str(self.search_term),
                sequence_type=sequence_operators.SequenceType.PROTEIN)

        if self.query_type == "seqmotif":
            return seqmotif_operators.SeqMotifOperator(
                pattern=str(self.search_term),
                sequence_type=seqmotif_operators.SequenceType.PROTEIN,
                pattern_type=seqmotif_operators.PatternType.SIMPLE)

        if self.query_type == "structure":
            return structure_operators.StructureOperator(
                pdb_entry_id=str(self.search_term),
                assembly_id=1,
                search_mode=structure_operators.StructureSearchMode.
                RELAXED_SHAPE_MATCH)

        if self.query_type == "chemical":
            return chemical_operators.ChemicalOperator(
                descriptor=str(self.search_term))

        operator_info = self._SUBTYPE_ATTRIBUTES.get(self.query_subtype)
        if operator_info is None:
            # A bare "text" query with no subtype has no attribute to search
            # against, so fall back to a full-text search.
            return text_operators.DefaultOperator(value=str(self.search_term))

        operator_name, attribute = operator_info
        if operator_name == "in":
            # The `in` operator matches against any of a list of values, so a
            # search term that is already a list is passed through as-is.
            values = (list(self.search_term) if isinstance(
                self.search_term, (list, tuple)) else [self.search_term])
            return text_operators.InOperator(attribute=attribute,
                                             values=values)
        if operator_name == "contains_words":
            return text_operators.ContainsWordsOperator(
                attribute=attribute, value=str(self.search_term))
        return text_operators.ExactMatchOperator(attribute=attribute,
                                                 value=str(self.search_term))

    def _custom_scan_params_search(self):
        """Searches using an explicitly-provided `scan_params` JSON query."""
        query_text = json.dumps(self.scan_params, indent=4)
        response = http_requests.request_limited(
            self.url,
            rtype="POST",
            headers={"Content-Type": "application/json"},
            data=query_text)

        if response is None or response.status_code != 200:
            warnings.warn("Retrieval failed, returning None")
            return None

        response_val = json.loads(response.text)

        if self.return_type == "entry":
            return walk_nested_dict(response_val,
                                    "identifier",
                                    maxdepth=25,
                                    outputs=[])
        return response_val

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

        Returns
        -------

        A list of PDB IDs matching the query (or, for the `polymer_entity`
        return type, the raw JSON response). Returns None if the request fails.
        """
        if self._uses_custom_scan_params:
            return self._custom_scan_params_search()

        try:
            return search_client.perform_search(
                search_operator=self._build_search_operator(),
                return_type=self._get_return_type_enum(),
                return_raw_json_dict=(self.return_type != "entry"),
                verbosity=False)
        except requests.RequestException:
            warnings.warn("Retrieval failed, returning None")
            return None


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
        return None

    return pdb_client.get_pdb_file(pdb_id, filetype_enum, compression)


# https://data.rcsb.org/migration-guide.html#chem-comp-description
def describe_chemical(chem_id):
    '''Look up the chemical description of a ligand

    Parameters
    ----------

    chem_id : string
        A 3 character string representing the full chemical sequence of
        interest (ie, NAG)

    Returns
    -------

    out : dict
        A dictionary containing the chemical description associated with
        the PDB ID

    Examples
    --------
    >>> chem_desc = describe_chemical('NAG')
    >>> print(chem_desc["rcsb_chem_comp_descriptor"]["smiles"])
    'CC(=O)NC1C(C(C(OC1O)CO)O)O'

    '''
    if (len(chem_id) > 3):
        raise Exception("Ligand id with more than 3 characters provided")

    return get_info(chem_id, url_root = 'https://data.rcsb.org/rest/v1/core/chemcomp/')

def get_ligands(pdb_id):
    """Return ligands of given PDB ID (DEPRECATED)"""
    warnings.warn(
        "get_ligands() is deprecated. Use pypdb.clients.data.DataFetcher with DataType.ENTRY "
        "and query for nonpolymer_entity_instances or related chemical component data.",
        DeprecationWarning,
        stacklevel=2
    )
    return None

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


def count_results(search_term, query_type="full_text"):
    '''Count the entries matching a search, without fetching them

    This is much faster than running the search and taking the length of the
    result, because RCSB returns only the tally.

    Parameters
    ----------

    search_term : str
        The search keyword, interpreted according to `query_type`

    query_type : str
        The type of query to perform (see help(Query) for the options)

    Returns
    -------

    out : int
        The number of entries matching the search

    Examples
    --------

    >>> print(count_results('ribosome'))
    9560

    '''
    query = Query(search_term, query_type=query_type)
    return search_client.perform_search(
        search_operator=query._build_search_operator(),
        return_type=query._get_return_type_enum(),
        request_options=search_client.RequestOptions(sort_by=None,
                                                     desc=None,
                                                     return_counts=True),
        verbosity=False)


def get_top_results(search_term, max_results=10, query_type="full_text"):
    '''Return the highest-scoring results of a search

    Fetching only the top hits is much faster than requesting every match.

    Parameters
    ----------

    search_term : str
        The search keyword, interpreted according to `query_type`

    max_results : int
        The maximum number of results to return

    query_type : str
        The type of query to perform (see help(Query) for the options)

    Returns
    -------

    out : list of str
        A descending-score list of the PDB IDs matching the search

    Examples
    --------

    >>> print(get_top_results('crispr', max_results=3))
    ['5FCL', '4XTK', '2Y9H']

    '''
    query = Query(search_term, query_type=query_type)
    return search_client.perform_search(
        search_operator=query._build_search_operator(),
        return_type=query._get_return_type_enum(),
        request_options=search_client.RequestOptions(
            result_start_index=0, num_results=max_results),
        verbosity=False)


def find_ligands(chem_id, max_results=None):
    '''Return the PDB entries that contain a given ligand

    Parameters
    ----------

    chem_id : str
        The chemical component ID of the ligand (e.g. 'NAG')

    max_results : int
        The maximum number of results to return. Defaults to all of them.

    Returns
    -------

    out : list of str
        A list of PDB IDs of entries containing the ligand

    Examples
    --------

    >>> print(find_ligands('NAG', max_results=3))
    ['10FT', '10GH', '10IC']

    '''
    if max_results is None:
        return Query(chem_id, query_type="chemical").search()
    return get_top_results(chem_id,
                           max_results=max_results,
                           query_type="chemical")


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
