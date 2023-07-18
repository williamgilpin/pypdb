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
- PubMed integrated data
- UniProt integrated data
- DrugBank integrated data
"""
from pypdb.clients.data.graphql.graphql import search_graphql

# TODO: convert to dataclass
class DataType:
    """
    General class that will host various data types, as detailed above.
    """
    def __init__(self, id):
        self.id = id

        self.properties = None
        self.json_query = None
        self.response = None

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

        # init self.properties to empty dict if None
        if self.properties is None:
            self.properties = {}

        # add properties to the dict
        for key, value in property.items():
            if key not in self.properties:
                self.properties[key] = value
            else:
                self.properties[key] += value
                self.properties[key] = list(set(self.properties[key]))

    def fetch_data(self):
        """
        Once the JSON query is created, fetch data from the PDB, using graphql.
        """
        if self.json_query is None:
            self.generate_json_query()

        response = search_graphql(self.json_query)
        self.response = response

class Entry(DataType):
    """
    DataType for Entry.

    https://data.rcsb.org/rest/v1/schema/entry
    """

    def check_pdb_id(self):
        """
        Check to see if we have valid pdb ids.
        """
        if isinstance(self.id, str):
            if len(self.id) != 4:
                raise ValueError
            self.id = [self.id]
        elif isinstance(self.id, list):
            if not all([isinstance(pid, str) for pid in self.id]):
                raise TypeError
            if not all([len(pid)==4 for pid in self.id]):
                raise ValueError
        else:
            raise TypeError

    def generate_json_query(self):
        """
        Given pdb id(s), and properties to fetch, generate json query.
        """
        if self.properties is None:
            print(f"ERROR: no properties given to generate JSON query.")
            raise ValueError

        self.check_pdb_id()

        if len(self.id) == 1:
            # we have a single pdb_id
            entry_str = f"entry(entry_id: \"{self.id[0]}\")"
        else:
            # we have a list of pdb_ids
            entry_str = "entries(entry_ids: [" + ",".join(f"\"{w}\"" for w in self.id) + "])"

        props_string = ""
        for key, val in self.properties.items():
            props_string += f"{key} {{" + ",".join(val) + "}"

        self.json_query = {'query': "{" + entry_str + "{" + props_string + "}}"}