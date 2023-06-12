from pypdb.clients.fasta.fasta_client import get_fasta_from_rcsb_entry
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators.sequence_operators import SequenceOperator

# Fetches the first sequence in the "6TML" fasta file
fasta_sequence = get_fasta_from_rcsb_entry("6TML", verbosity=False)[0].sequence

# Performs sequence search ('BLAST'-like) using the FASTA sequence
results = perform_search(
        return_type=ReturnType.ENTRY,
        verbosity=False,
    search_operator=SequenceOperator(
        sequence=fasta_sequence,
        identity_cutoff=0.99,
        evalue_cutoff=1000
      ),
    return_with_scores=True
)