# PyPDB Text Search

## Helpful Links

The Search logic here is a Python wrapper around the RCSB's search logic.
For in-the-weeds details on how each operator works, prefer to look at the
[RCSB Search API documentation](https://search.rcsb.org/index.html)

The search operators defined within the `operators` directory support querying
RCSB attributes against the appropriate `SearchService`. For example, if
you are querying the RCSB Text Search Service (`SearchService.TEXT`), all
operators within `text_operators.py` should be supported.

For a list of RCSB attributes associated with structures you can search, see
[RCSB's List of Attributes to Search](http://search.rcsb.org/search-attributes.html).
Note that not every structure will have every attribute.

Two querying functions are currently supported by PyPDB:

* `perform_search`: This function is good for simple queries
* `perform_search_with_graph`: This function allows building complicated queries using RCSB's query node syntax.

## `perform_search` Examples

### Search for all entries that mention the word 'ribosome'
```
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators

search_service = SearchService.TEXT
search_operator = text_operators.DefaultOperator(value="ribosome")
return_type = ReturnType.ENTRY

results = perform_search(search_service, search_operator, return_type)
```

### Search for polymers from 'Mus musculus'
```
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators

search_service = SearchService.TEXT
search_operator = text_operators.ExactMatchOperator(value="Mus musculus",
                                                    attribute="rcsb_entity_source_organism.taxonomy_lineage.name")
return_type = ReturnType.POLYMER_ENTITY

results = perform_search(search_service, search_operator, return_type)
```

### Search for non-polymers from 'Mus musculus' or 'Homo sapiens'
```
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators

search_operator = text_operators.InOperator(values=["Mus musculus", "Homo sapiens"],
                                            attribute="rcsb_entity_source_organism.taxonomy_lineage.name")
return_type = ReturnType.NON_POLYMER_ENTITY

results = perform_search(search_service, search_operator, return_type)
```

### Search for polymer instances whose titles contain "actin" or "binding" or "protein"
```
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators

search_service = SearchService.TEXT
search_operator = text_operators.ContainsWordsOperator(value="actin-binding protein",
                                            attribute="struct.title")
return_type = ReturnType.POLYMER_INSTANCE

results = perform_search(search_service, search_operator, return_type)
```

### Search for assemblies that contain the words "actin binding protein"
(must be in that order).

For example, "actin-binding protein" and "actin binding protein" will match,
but "protein binding actin" will not.
```
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators

search_service = SearchService.TEXT
search_operator = text_operators.ContainsPhraseOperator(value="actin-binding protein",
                                            attribute="struct.title")
return_type = ReturnType.ASSEMBLY

results = perform_search(search_service, search_operator, return_type)
```

### Search for entries released in 2019 or later
```
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators

search_service = SearchService.TEXT
search_operator = text_operators.ComparisonOperator(
       value="2019-01-01T00:00:00Z",
       attribute="rcsb_accession_info.initial_release_date",
       comparison_type=text_operators.ComparisonType.GREATER)
return_type = ReturnType.ENTRY

results = perform_search(search_service, search_operator, return_type)
```

### Search for entries released only in 2019 or later
```
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators

search_service = SearchService.TEXT
search_operator = text_operators.RangeOperator(
    from_value="2019-01-01T00:00:00Z",
    to_value="2020-01-01T00:00:00Z",
    include_lower=True,
    include_upper=False,
    attribute="rcsb_accession_info.initial_release_date")
return_type = ReturnType.ENTRY

results = perform_search(search_service, search_operator, return_type)
```
### Search for structures under 4 angstroms of resolution
```
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators

search_service = SearchService.TEXT
search_operator = text_operators.ComparisonOperator(
           value=4,
           attribute="rcsb_entry_info.resolution_combined",
           comparison_type=text_operators.ComparisonType.LESS)
return_type = ReturnType.ENTRY

results = perform_search(search_service, search_operator, return_type)
```


### Search for structures with a given attribute.

(Admittedly every structure has a release date, but the same logic would
 apply for a more sparse RCSB attribute).

```
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators

search_service = SearchService.TEXT
search_operator = text_operators.ExistsOperator(
    attribute="rcsb_accession_info.initial_release_date")
return_type = ReturnType.ENTRY

results = perform_search(search_service, search_operator, return_type)
```

### Search for top 100 structures matching the given protein sequence, by date

(this sequence matches the SARS-CoV-2 NSP3 macrodomain)

```
from pypdb.clients.search.search_client import perform_search, RequestOptions
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators.sequence_operators import SequenceOperator
from pypdb.clients.search.operators.sequence_operators import SequenceType

results = perform_search(
    search_service=SearchService.SEQUENCE,
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
```
from pypdb.clients.fasta.fasta_client import get_fasta_from_rcsb_entry
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators.sequence_operators import SequenceOperator

# Fetches the first sequence in the "6TML" fasta file
fasta_sequence = get_fasta_from_rcsb_entry("6TML")[0].sequence

# Performs sequence search ('BLAST'-like) using the FASTA sequence
results = perform_search(
    search_service=SearchService.SEQUENCE,
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

```
from pypdb.clients.search.search_client import perform_search_with_graph
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.search_client import QueryNode, QueryGroup, LogicalOperator
from pypdb.clients.search.operators import text_operators

# QueryNode associated with structures with under 4 Angstroms of resolution
under_4A_resolution_operator = text_operators.ComparisonOperator(
       value=4,
       attribute="rcsb_entry_info.resolution_combined",
       comparison_type=text_operators.ComparisonType.GREATER)
under_4A_query_node = QueryNode(SearchService.TEXT,
                                  under_4A_resolution_operator)

# QueryNode associated with entities containing 'Mus musculus' lineage
is_mus_operator = text_operators.ExactMatchOperator(
            value="Mus musculus",
            attribute="rcsb_entity_source_organism.taxonomy_lineage.name")
is_mus_query_node = QueryNode(SearchService.TEXT, is_mus_operator)

# QueryNode associated with entities containing 'Homo sapiens' lineage
is_human_operator = text_operators.ExactMatchOperator(
            value="Homo sapiens",
            attribute="rcsb_entity_source_organism.taxonomy_lineage.name")
is_human_query_node = QueryNode(SearchService.TEXT, is_human_operator)

# QueryGroup associated with being either human or `Mus musculus`
is_human_or_mus_group = QueryGroup(
    queries = [is_mus_query_node, is_human_query_node],
    logical_operator = LogicalOperator.OR
)

# QueryGroup associated with being ((Human OR Mus) AND (Under 4 Angstroms))
is_under_4A_and_human_or_mus_group = QueryGroup(
    queries = [is_human_or_mus_group, under_4A_query_node],
    logical_operator = LogicalOperator.AND
)

return_type = ReturnType.ENTRY

results = perform_search_with_graph(
  query_object=is_under_4A_and_human_or_mus_group,
  return_type=return_type)
```

## Search for Calcium-Bound Calmodulin Structures

Note that "1CLL" corresponds to a Calmodulin structure bound to Ca2+.

Also, searching for `rcsb_chem_comp_container_identifiers.comp_id` with
an exact match to `"CA"` yields only structures in complex with Ca2+
(filtering out structures in complex with other metals like strontium).

```
from pypdb.clients.search.search_client import perform_search_with_graph
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.search_client import QueryNode, QueryGroup, LogicalOperator
from pypdb.clients.search.operators import text_operators, structure_operators

is_similar_to_1CLL = QueryNode(
  search_service=SearchService.STRUCTURE,
  search_operator=structure_operators.StructureOperator(
      pdb_entry_id="1CLL",
      assembly_id=1,
      search_mode=structure_operators.StructureSearchMode.STRICT_SHAPE_MATCH
  )
)

is_in_complex_with_calcium = QueryNode(
  search_service=SearchService.TEXT,
  search_operator=text_operators.ExactMatchOperator(
    attribute="rcsb_chem_comp_container_identifiers.comp_id",
    value="CA"
  )
)

results = perform_search_with_graph(
  query_object=QueryGroup(
    logical_operator=LogicalOperator.AND,
    queries=[is_similar_to_1CLL, is_in_complex_with_calcium]
  ),
  return_type=ReturnType.ENTRY
)
