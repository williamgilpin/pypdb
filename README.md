A simple Python API for submitting and processing queries to the RCSB Protein Data Bank (PDB)

This program uses the requests library as well as several xml parsing libraries in order to create calls to the PDB database, which has a REST API that plays well with xml requests. The returned data has an elaborate structure due to the format of the query responses; future versions will work on representing this information with some sort of high-level Python object.



William Gilpin, 2015