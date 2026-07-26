# PyPDB Text Search

## Helpful Links

The Search logic here is a Python wrapper around the RCSB's search logic.
For in-the-weeds details on how each operator works, prefer to look at the
[RCSB Search API documentation](https://search.rcsb.org/index.html)

The search operators defined within the `operators` directory support querying
RCSB attributes against the appropriate `if
you are querying the RCSB Text Search Service (`all
operators within `text_operators.py` should be supported.

For a list of RCSB attributes associated with structures you can search, see
[RCSB's List of Structure Attributes to Search](https://search.rcsb.org/structure-search-attributes.html) and [RCSB's List of Chemical Attributes to Search](https://search.rcsb.org/chemical-search-attributes.html)
Note that not every structure will have every attribute.

Two querying functions are currently supported by PyPDB:

* `perform_search`: This function is good for simple queries
* `perform_search_with_graph`: This function allows building complicated queries using RCSB's query node syntax.

## `perform_search` Examples

### Search for all entries that mention the word 'ribosome'

```python
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators


search_operator = text_operators.DefaultOperator(value="ribosome")
return_type = ReturnType.ENTRY

results = perform_searchsearch_operator, return_type)
```

### Search for polymers from 'Mus musculus'

```python
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators


search_operator = text_operators.ExactMatchOperator(value="Mus musculus",
                                                    attribute="rcsb_entity_source_organism.taxonomy_lineage.name")
return_type = ReturnType.POLYMER_ENTITY

results = perform_search(search_operator, return_type)
```

### Search for non-polymers from 'Mus musculus' or 'Homo sapiens'

```python
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators

search_operator = text_operators.InOperator(values=["Mus musculus", "Homo sapiens"],
                                            attribute="rcsb_entity_source_organism.taxonomy_lineage.name")
return_type = ReturnType.NON_POLYMER_ENTITY

results = perform_search(search_operator, return_type)
```

### Search for polymer instances whose titles contain "actin" or "binding" or "protein"

```python
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators


search_operator = text_operators.ContainsWordsOperator(value="actin-binding protein",
                                            attribute="struct.title")
return_type = ReturnType.POLYMER_INSTANCE

results = perform_search(search_operator, return_type)
```

### Search for assemblies that contain the words "actin binding protein"

(must be in that order).

For example, "actin-binding protein" and "actin binding protein" will match,
but "protein binding actin" will not.

```python
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators


search_operator = text_operators.ContainsPhraseOperator(value="actin-binding protein",
                                            attribute="struct.title")
return_type = ReturnType.ASSEMBLY

results = perform_search(search_operator, return_type)
```

### Search for entries released in 2019 or later

```python
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators


search_operator = text_operators.ComparisonOperator(
       value="2019-01-01T00:00:00Z",
       attribute="rcsb_accession_info.initial_release_date",
       comparison_type=text_operators.ComparisonType.GREATER)
return_type = ReturnType.ENTRY

results = perform_search(search_operator, return_type)
```

### Search for entries released only in 2019 or later

```python
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators


search_operator = text_operators.RangeOperator(
    from_value="2019-01-01T00:00:00Z",
    to_value="2020-01-01T00:00:00Z",
    include_lower=True,
    include_upper=False,
    attribute="rcsb_accession_info.initial_release_date")
return_type = ReturnType.ENTRY

results = perform_search(search_operator, return_type)
```

### Search for structures under 4 angstroms of resolution

```python
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators


search_operator = text_operators.ComparisonOperator(
           value=4,
           attribute="rcsb_entry_info.resolution_combined",
           comparison_type=text_operators.ComparisonType.LESS)
return_type = ReturnType.ENTRY

results = perform_search(search_operator, return_type)
```

### Search for structures with a given attribute

(Admittedly every structure has a release date, but the same logic would
 apply for a more sparse RCSB attribute).

```python
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators


search_operator = text_operators.ExistsOperator(
    attribute="rcsb_accession_info.initial_release_date")
return_type = ReturnType.ENTRY

results = perform_search(search_operator, return_type)
```

### Search for top 100 structures matching the given protein sequence, by date

(this sequence matches the SARS-CoV-2 NSP3 macrodomain)

```python
from pypdb.clients.search.search_client import perform_search, RequestOptions
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators.sequence_operators import SequenceOperator
from pypdb.clients.search.operators.sequence_operators import SequenceType

results = perform_search(
        return_type=ReturnType.ENTRY,
    search_operator=SequenceOperator(
        sequence_type=SequenceType.PROTEIN, # if not explicitly specified, this will autoresolve
        sequence=(
          "SMVNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDY"
          "IATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADP"
          "IHSLRVCVDTVRTNVYLAVFDKNLYDKLVSSFL"),
        identity_cutoff=0.99,
        evalue_cutoff=1000
      ),
    request_options=RequestOptions(
        result_start_index=0,
        num_results=100,
        sort_by="rcsb_accession_info.initial_release_date",
        desc=False
      ),
    return_with_scores=True
)
```

### Search for structures that match the sequence of an existing RCSB entry

```python
from pypdb.clients.fasta.fasta_client import get_fasta_from_rcsb_entry
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators.sequence_operators import SequenceOperator

# Fetches the first sequence in the "6TML" fasta file
fasta_sequence = get_fasta_from_rcsb_entry("6TML", verbosity=True)[0].sequence

# Performs sequence search ('BLAST'-like) using the FASTA sequence
results = perform_search(
        return_type=ReturnType.ENTRY,
    search_operator=SequenceOperator(
        sequence=fasta_sequence,
        identity_cutoff=0.99,
        evalue_cutoff=1000
      ),
    return_with_scores=True
)
```

## `perform_search_with_graph` Example

### Search for 'Mus musculus' or 'Homo sapiens' structures after 2019

```python
from pypdb.clients.search.search_client import perform_search_with_graph
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.search_client import QueryGroup, LogicalOperator
from pypdb.clients.search.operators import text_operators

# SearchOperator associated with structures with under 4 Angstroms of resolution
under_4A_resolution_operator = text_operators.ComparisonOperator(
       value=4,
       attribute="rcsb_entry_info.resolution_combined",
       comparison_type=text_operators.ComparisonType.GREATER)

# SearchOperator associated with entities containing 'Mus musculus' lineage
is_mus_operator = text_operators.ExactMatchOperator(
            value="Mus musculus",
            attribute="rcsb_entity_source_organism.taxonomy_lineage.name")

# SearchOperator associated with entities containing 'Homo sapiens' lineage
is_human_operator = text_operators.ExactMatchOperator(
            value="Homo sapiens",
            attribute="rcsb_entity_source_organism.taxonomy_lineage.name")

# QueryGroup associated with being either human or `Mus musculus`
is_human_or_mus_group = QueryGroup(
    queries = [is_mus_operator, is_human_operator],
    logical_operator = LogicalOperator.OR
)

# QueryGroup associated with being ((Human OR Mus) AND (Under 4 Angstroms))
is_under_4A_and_human_or_mus_group = QueryGroup(
    queries = [is_human_or_mus_group, under_4A_resolution_operator],
    logical_operator = LogicalOperator.AND
)

results = perform_search_with_graph(
  query_object=is_under_4A_and_human_or_mus_group,
  return_type=ReturnType.ENTRY)
```

## Search for Calcium-Bound Calmodulin Structures

Note that "1CLL" corresponds to a Calmodulin structure bound to Ca2+.

Also, searching for `rcsb_chem_comp_container_identifiers.comp_id` with
an exact match to `"CA"` yields only structures in complex with Ca2+
(filtering out structures in complex with other metals like strontium).

```python
from pypdb.clients.search.search_client import perform_search_with_graph
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.search_client import QueryGroup, LogicalOperator
from pypdb.clients.search.operators import text_operators, structure_operators

is_similar_to_1CLL = structure_operators.StructureOperator(
    pdb_entry_id="1CLL",
    assembly_id=1,
    search_mode=structure_operators.StructureSearchMode.STRICT_SHAPE_MATCH
)

is_in_complex_with_calcium = text_operators.ExactMatchOperator(
    attribute="rcsb_chem_comp_container_identifiers.comp_id",
    value="CA"
)

results = perform_search_with_graph(
  query_object=QueryGroup(
    logical_operator=LogicalOperator.AND,
    queries=[is_similar_to_1CLL, is_in_complex_with_calcium]
  ),
  return_type=ReturnType.ENTRY
)

## Count results without fetching them

Setting `return_counts` returns only the number of matching entries, which is
much faster than retrieving every hit.

```python
from pypdb.clients.search.search_client import perform_search, ReturnType, RequestOptions
from pypdb.clients.search.operators import text_operators

count = perform_search(
    search_operator=text_operators.DefaultOperator(value="ribosome"),
    return_type=ReturnType.ENTRY,
    request_options=RequestOptions(sort_by=None, desc=None, return_counts=True)
)
print(count)  # e.g. 9560
```

## Fetch only the top results

```python
from pypdb.clients.search.search_client import perform_search, ReturnType, RequestOptions
from pypdb.clients.search.operators import text_operators

results = perform_search(
    search_operator=text_operators.DefaultOperator(value="crispr"),
    return_type=ReturnType.ENTRY,
    request_options=RequestOptions(result_start_index=0, num_results=10)
)
```

## Include computed structure models

By default only experimentally-determined structures are searched. Computed
models (such as AlphaFold predictions) can be included as well.

```python
from pypdb.clients.search.search_client import (
    perform_search, ReturnType, RequestOptions, ResultsContentType)
from pypdb.clients.search.operators import text_operators

results = perform_search(
    search_operator=text_operators.DefaultOperator(value="ribosome"),
    return_type=ReturnType.ENTRY,
    request_options=RequestOptions(
        result_start_index=0,
        num_results=10,
        results_content_type=[ResultsContentType.EXPERIMENTAL,
                              ResultsContentType.COMPUTATIONAL])
)
print(results)  # includes IDs like 'AF_AFA0A0K0DUR6F1'
```

## Search molecular definitions by chemical formula

The `mol_definition` return type returns chemical component (and BIRD) IDs,
rather than structures.

```python
from pypdb.clients.search.search_client import perform_search, ReturnType
from pypdb.clients.search.operators.chemical_operators import ChemicalFormulaOperator

results = perform_search(
    search_operator=ChemicalFormulaOperator(formula="C12H17N4OS"),
    return_type=ReturnType.MOL_DEFINITION
)
print(results)  # ['VIB'] (thiamine)
```

## Search molecular definitions by attribute

Attribute searches against chemical definitions use the operators in
`selection_operators`, which are dispatched to RCSB's `text_chem` service.

```python
from pypdb.clients.search.search_client import perform_search, ReturnType
from pypdb.clients.search.operators import selection_operators

results = perform_search(
    search_operator=selection_operators.ChemicalExactMatchOperator(
        attribute="rcsb_chem_comp_container_identifiers.comp_id",
        value="NAG"),
    return_type=ReturnType.MOL_DEFINITION
)
print(results)  # ['NAG']
```

## Search by Enzyme Classification (EC) number

EC numbers are stored as a lineage, so a partial EC number matches every class
beneath it in the hierarchy (`"2.3.2"` returns all of its sub-classes).

```python
from pypdb import Query

found_pdbs = Query("2.3.2.5", query_type="ec_number").search()
```

The equivalent search using the operator API, which also accepts several EC
numbers at once:

```python
from pypdb.clients.search.search_client import perform_search, ReturnType
from pypdb.clients.search.operators import text_operators

results = perform_search(
    search_operator=text_operators.InOperator(
        attribute="rcsb_polymer_entity.rcsb_ec_lineage.id",
        values=["2.3.2.5"]),
    return_type=ReturnType.ENTRY
)
```

## Search for a 3D structure motif

A structure motif search finds spatial arrangements of residues, such as a
catalytic triad, regardless of the sequence or fold they appear in. The motif
is defined by pointing at between 2 and 10 residues of an existing entry,
addressed by their mmCIF `label_asym_id` and `label_seq_id`.

```python
from pypdb.clients.search.search_client import perform_search, ReturnType
from pypdb.clients.search.operators.strucmotif_operators import (
    AtomPairingScheme, StructureMotifResidue, StrucMotifOperator)

# The enolase superfamily's catalytic residues, as they sit in 2MNR
results = perform_search(
    search_operator=StrucMotifOperator(
        pdb_entry_id="2MNR",
        residues=[
            StructureMotifResidue(label_asym_id="A", label_seq_id=162),
            StructureMotifResidue(label_asym_id="A", label_seq_id=193),
            StructureMotifResidue(label_asym_id="A", label_seq_id=219),
        ],
        rmsd_cutoff=2,
        atom_pairing_scheme=AtomPairingScheme.ALL),
    return_type=ReturnType.ENTRY
)
```

Individual positions can permit other residue types via `exchanges`, which
broadens the search:

```python
StructureMotifResidue(label_asym_id="A", label_seq_id=162,
                      exchanges=["LYS", "HIS"])
```

## Search for chemically similar molecules

`ChemicalOperator` matches molecular graphs exactly or as substructures.
To rank molecular definitions by chemical similarity instead, use
`ChemicalSimilarityOperator`:

```python
from pypdb.clients.search.search_client import perform_search, ReturnType
from pypdb.clients.search.operators.chemical_operators import ChemicalSimilarityOperator

results = perform_search(
    search_operator=ChemicalSimilarityOperator(
        descriptor="CC(=O)NC1C(O)OC(CO)C(O)C1O"),
    return_type=ReturnType.MOL_DEFINITION
)
```

## Tally results with facets

A facet buckets the matching results by an attribute's value, and is returned
alongside the hits, so counting how many structures used each experimental
method takes one request rather than one per method.

```python
from pypdb.clients.search.search_client import (
    Facet, RequestOptions, ReturnType, perform_search)
from pypdb.clients.search.operators import text_operators

response = perform_search(
    search_operator=text_operators.DefaultOperator(value="hemoglobin"),
    return_type=ReturnType.ENTRY,
    request_options=RequestOptions(
        result_start_index=0,
        num_results=1,
        facets=[Facet(name="Methods", attribute="exptl.method")]),
    # Facets are only present in the raw response, not in the ID list
    return_raw_json_dict=True
)

for bucket in response["facets"][0]["buckets"]:
    print(bucket["label"], bucket["population"])
# X-RAY DIFFRACTION 8253
# ELECTRON MICROSCOPY 643
# SOLUTION NMR 151
```

## Group redundant results

Grouping collapses related hits together — for example, one group per unique
sequence instead of every redundant copy of it.

```python
from pypdb.clients.search.search_client import (
    GroupBy, GroupByReturnType, RequestOptions, ReturnType, perform_search)
from pypdb.clients.search.operators import text_operators

response = perform_search(
    search_operator=text_operators.DefaultOperator(value="hemoglobin"),
    return_type=ReturnType.POLYMER_ENTITY,
    request_options=RequestOptions(
        result_start_index=0,
        num_results=1,
        group_by=GroupBy(aggregation_method="sequence_identity",
                         similarity_cutoff=100),
        group_by_return_type=GroupByReturnType.GROUPS),
    return_raw_json_dict=True
)

print(response["group_by_count"])  # number of distinct groups
```

Use `GroupByReturnType.REPRESENTATIVES` to get a single member of each group
instead of the whole group.
