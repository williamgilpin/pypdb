# PyPDB Data Fetch from the PDB

## Helpful Links

The data fetch module here is a Python wrapper for the [graphQL API](https://data.rcsb.org/#fetch-data-graphql).
PDB's data API [organizes the data in the following way](https://data.rcsb.org/#data-organization):

* `entry`
* `entity`
    * `polymer_entity`
    * `branched_entity`
    * `nonpolymer_entity`
* `entity_instance`
    * `polymer_entity_instance`
    * `branched_entity_instance`
    * `nonpolymer_entity_instance`
* `assembly`
* `chemical_component`

In addition to these, the following are also options in the PDB, but are currently not implemented in PyPDB:

* `PubMed`
* `UniProt`
* `DrugBank`

The data schemas for all of these data types can be viewed [here](https://data.rcsb.org/#data-schema).
These schemas allow the user to determine what keywords to ask for.
The queries can be tested in-browser using the [GraphiQL tool](https://data.rcsb.org/graphql/index.html?query=%7B%0A%20%20entries(entry_ids%3A%20%5B%224HHB%22%5D)%20%7B%0A%20%20%20%20rcsb_id%0A%20%20%20%20struct%20%7B%0A%20%20%20%20%20%20title%0A%20%20%20%20%7D%0A%20%20%20%20exptl%20%7B%0A%20%20%20%20%20%20method%0A%20%20%20%20%7D%0A%20%20%7D%0A%7D).

## Examples

All of the functionaly, and thus examples below, require the following imports:

```python
from pypdb.clients.data.data_types import DataFetcher, DataType
```

### Fetch entries using PDB IDs

If we want to fetch some information about the PDB entries `4HHB`, `12CA`, and `3PQR`, we first create an instance of `DataFetcher`:

```python
entry = DataFetcher(["4HHB", "12CA", "3PQR"], DataType.ENTRY)
```

The properties we will fetch for needs to be given as a python dictionary, commensurate with the [data schemas](https://data.rcsb.org/#data-schema):

```python
property = {"exptl": ["method", "details"], "cell":["length_a", "length_b", "length_c"]}
entry.add_property(property)
```

Then we fetch the data

```python
entry.fetch_data()
```

where `entry.response` now contains a Python dictionary generated from the JSON formatted information fetched from the PDB.
It is possible to convert this to a Pandas dataframe:

```python
df = entry.return_data_as_pandas_df()
```

### Fetch Assemblies

Similarly to the `entry` case:

```python
assembly = DataFetcher(["4HHB-1", "12CA-1", "3PQR-1"], DataType.ASSEMBLY)
property = {"rcsb_assembly_info": ["entry_id", "assembly_id", "polymer_entity_instance_count"]}

assembly.add_property(property)
assembly.fetch_data()
```

Note that the IDs provided must be of the form `[entry_id]-[assembly_id]`. 

### Fetch Polymer Entities

```python
fetcher = DataFetcher(["2CPK_1","3WHM_1","2D5Z_1"], DataType.POLYMER_ENTITY)
property = {"rcsb_id": [], 
            "rcsb_entity_source_organism": ["ncbi_taxonomy_id", "ncbi_scientific_name"],
            "rcsb_cluster_membership": ["cluster_id", "identity"]}

fetcher.add_property(property)
fetcher.fetch_data()
```

The IDs provided must be of the form `[entry_id]_[entity_id]`.

### Fetch Polymer Entity Instance

```python
fetcher = DataFetcher(["4HHB.A", "12CA.A", "3PQR.A"], DataType.POLYMER_ENTITY_INSTANCE)
property = {"rcsb_id": [],
            "rcsb_polymer_instance_annotation": ["annotation_id", "name", "type"]}
fetcher.add_property(property)
fetcher.fetch_data()
```

In this case, IDs are of the form `[entry_id].[entity_id]`.

### Fetch Branched Entity

```python
fetcher = DataFetcher(["5FMB_2", "6L63_3"], DataType.BRANCHED_ENTITY)
property = {"pdbx_entity_branch": ["type"],
            "pdbx_entity_branch_descriptor": ["type", "descriptor"]}

fetcher.add_property(property)
fetcher.fetch_data()
```

### Fetch Chemical Components

```python
fetcher = DataFetcher(["NAG","EBW"], DataType.CHEMICAL_COMPONENT)
property = {"rcsb_id":[], "chem_comp": ["type", "formula_weight","name","formula"],
            "rcsb_chem_comp_info":["initial_release_date"]}
fetcher.add_property(property)
fetcher.fetch_data()
```