"""
Class for data types that can be accessed in the PDB DATA-API
https://data.rcsb.org/#data-organization

Namely:
- entry
- polymer entity
- branched entity
- non-polymer entity
- polymer instance
- branched instance
- non-polymer instance
- assembly
- chemical component
(currently not implemented:)
- PubMed integrated data
- UniProt integrated data
- DrugBank integrated data
"""
from dataclasses import dataclass, field
from enum import Enum

#TODO: handle batch requests

from pypdb.clients.data.graphql.graphql import search_graphql

class DataType(Enum):
    ENTRY = "entries"
    POLYMER_ENTITY = "polymer_entities"
    BRANCHED_ENTITY = "branched_entities"
    NONPOLYMER_ENTITY = "nonpolymer_entities"
    POLYMER_ENTITY_INSTANCE = "polymer_entity_instances"
    BRANCHED_ENTITY_INSTANCE = "branched_entity_instances"
    NONPOLYMER_ENTITY_INSTANCE = "nonpolymer_entity_instances"
    ASSEMBLY = "assemblies"
    CHEMICAL_COMPONENT = "chem_comps"
    INTERFACE = "interfaces"
    # PubMed and UniProt integrated data are only queryable one ID at a time
    PUBMED = "pubmed"
    UNIPROT = "uniprot"


# Data types that the PDB only exposes as single-ID (rather than batch) queries
SINGULAR_DATA_TYPES = {DataType.PUBMED, DataType.UNIPROT}

# GraphQL argument name used to pass IDs, per data type
_ID_ARGUMENT_NAMES = {
    DataType.ENTRY: "entry_ids",
    DataType.ASSEMBLY: "assembly_ids",
    DataType.CHEMICAL_COMPONENT: "comp_ids",
    DataType.INTERFACE: "interface_ids",
    DataType.PUBMED: "pubmed_id",
    DataType.UNIPROT: "uniprot_id",
}

@dataclass
class DataFetcher:
    """
    General class that will host various data types, as detailed above.
    """

    id: str | int | list
    data_type: DataType

    properties: dict = field(default_factory=dict)
    json_query: dict = field(default_factory=dict)
    response: dict = field(default_factory=dict)

    def __post_init__(self):
        """
        Check types of IDs given, format accordingly.
        """

        # PubMed IDs are numeric, so accept ints alongside string identifiers
        if isinstance(self.id, (str, int)):
            self.id = [self.id]

        if self.data_type in SINGULAR_DATA_TYPES:
            if len(self.id) > 1:
                print(f"WARNING: {self.data_type.value} accepts only one ID "
                      "at a time; using the first one.")
                self.id = self.id[:1]
        elif "entit" in self.data_type.value and "instance" not in self.data_type.value:
            for id in self.id:
                if '_' not in id:
                    print(f"WARNING: {id} not valid for {self.data_type.value}.")
        elif "instance" in self.data_type.value:
            for id in self.id:
                if '.' not in id:
                    print(f"WARNING: {id} not valid for {self.data_type.value}.")
        elif self.data_type == DataType.ASSEMBLY:
            for id in self.id:
                if '-' not in id:
                    print(f"WARNING: {id} not valid for {self.data_type.value}.")

    def add_property(self, property):
        """
        Add property to the list of data to fetch from the PDB.

        property is a python dict, with keys as properties, and
        values as subproperties.

        e.g.:

        {"cell": ["volume", "angle_beta"], "exptl": ["method"]}

        If the user is trying to add a property that already exists,
        the subproperties are merged.
        """
        # check input data type
        if not isinstance(property, dict):
            raise TypeError
        # check data types of keys in dict
        if not all([isinstance(key, str) for key in property.keys()]):
            raise TypeError
        # check that values are lists of strings
        for key, value in property.items():
            if isinstance(value, str):
                property[key] = [value]
            elif not isinstance(value, list):
                raise TypeError
            else:
                if not all([isinstance(val, str) for val in value]):
                    raise TypeError

        # add properties to the dict
        for key, value in property.items():
            if key not in self.properties:
                self.properties[key] = value
            else:
                self.properties[key] += value
                self.properties[key] = list(set(self.properties[key]))

    def generate_json_query(self):
        """
        Given IDs, data type, and properties to fetch, create JSON query that
        will utilize graphql.
        """
        if not self.properties:
            print("ERROR: no properties given to generate JSON query.")
            raise ValueError

        q_str = _ID_ARGUMENT_NAMES.get(self.data_type)
        if q_str is None:
            if "instance" in self.data_type.value:
                q_str = "instance_ids"
            else:
                q_str = "entity_ids"

        if self.data_type in SINGULAR_DATA_TYPES:
            # Singular data types take a single (unquoted, for PubMed) ID
            # rather than a list.
            value = self.id[0]
            id_str = value if self.data_type == DataType.PUBMED else f"\"{value}\""
            data_str = f"{self.data_type.value}({q_str}: {id_str})"
        else:
            data_str = f"{self.data_type.value}({q_str}: [" + ",".join(
                f"\"{w}\"" for w in self.id) + "])"

        props_string = ""
        for key, val in self.properties.items():
            if len(val) == 0:
                props_string += f"{key},"
            else:
                props_string += f"{key} {{" + ",".join(val) + "}"

        self.json_query = {'query': "{" + data_str + "{" + props_string + "}}"}


    def fetch_data(self):
        """
        Once the JSON query is created, fetch data from the PDB, using graphql.
        """
        if not self.json_query:
            self.generate_json_query()

        response = search_graphql(self.json_query)

        if "errors" in response:
            print("ERROR encountered in fetch_data().")
            for error in response['errors']:
                print(error['message'])

            return

        self.response = response

        returned_data = self.response['data'][self.data_type.value]
        if self.data_type in SINGULAR_DATA_TYPES:
            if returned_data is None:
                print("WARNING: one or more IDs not found in the PDB.")
        elif len(returned_data) != len(self.id):
            print("WARNING: one or more IDs not found in the PDB.")

    def return_data_as_df_dict(self):
        """
        Return the fetched data as a dict usable by pandas or polars.
        """
        if not self.response:
            return None

        data = self.response['data'][self.data_type.value]

        if self.data_type in SINGULAR_DATA_TYPES:
            # Singular data types return one object rather than a list
            data = [data] if data is not None else []

        # flatten data dictionary by joining property and subproperty names
        data_flat = {}
        for i, entry in enumerate(data):
            id = self.id[i]
            curr_dict = {}
            for key, values in entry.items():
                if isinstance(values, list):
                    v = values[0]
                else:
                    v = values
                if isinstance(v, str):
                    new_key = f"{key}"
                    curr_dict[new_key] = v
                else:
                    for subprop, val in v.items():
                        new_key = f"{key}.{subprop}"
                        curr_dict[new_key] = val
            data_flat[id] = curr_dict

        return data_flat

    def return_data_as_pandas_df(self):
        """
        Backwards-compatible helper returning fetched data as a pandas DataFrame.
        """
        data = self.return_data_as_df_dict()
        if data is None:
            return None

        try:
            import pandas as pd
        except ImportError as exc:
            return data

        return pd.DataFrame.from_dict(data, orient="index")
