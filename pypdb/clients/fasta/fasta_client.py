"""Client to fetch FASTA files associated with structures from RCSB."""

from dataclasses import dataclass
import re
import requests
from typing import Dict, List
import warnings

FASTA_BASE_URL = "https://www.rcsb.org/fasta/entry/"

# Fasta Sequences are uniquely identified by a polymeric entity ID that looks
# like `${ENTRY_ID}_{SEQUENCE_NUMBER}` (e.g. `5JUP_1` or `6TML_10`)
PolymerEntity = str  # Defines type-alias (Polymer entity IDs are strings)


@dataclass
class FastaSequence:
    """Class containing data for one FASTA sequence (one of many in a file)."""
    # Polymeric entity ID uniquely identifying this sequence
    entity_id: PolymerEntity  # e.g. `"5RU3_1"`
    # Chains associated with this sequence
    chains: List[str]  # e.g. `["A", "B"]`
    # Sequence associated with this entity
    sequence: str
    # Un-processed FASTA header for a sequence
    # (e.g. `5RU3_1|Chains A,B|Non-structural protein 3|Severe acute respiratory syndrome coronavirus 2 (2697049)`)
    fasta_header: str


def _parse_fasta_text_to_list(raw_fasta_text: str) -> List[FastaSequence]:
    """Parses raw FASTA response into easy-to-use dict representation."""
    # Gets list of FASTA chunks (one per sequence)
    fasta_sequence_chunks = raw_fasta_text.strip().split(">")[1:]

    fasta_list = []
    for fasta_sequence_chunk in fasta_sequence_chunks:
        chunk_lines = fasta_sequence_chunk.split("\n")
        fasta_header = chunk_lines[0]
        fasta_sequence = "".join(chunk_lines[1:])

        header_segments = fasta_header.split("|")
        entity_id = header_segments[0]
        # Derives associated chains from header
        chains = re.sub("Chains? ", "", header_segments[1]).split(",")

        fasta_list.append(
            FastaSequence(entity_id=entity_id,
                          chains=chains,
                          sequence=fasta_sequence,
                          fasta_header=fasta_header))
    return fasta_list


def get_fasta_from_rcsb_entry(rcsb_id: str,
                              verbosity: bool = True,
                              ) -> List[FastaSequence]:
    """Fetches FASTA sequence associated with PDB structure from RCSB.

    Args:
      rcsb_id: RCSB accession code of the structure of interest. E.g. `"5RU3"`
      verbosity: Print out the search query to the console (default: True)

    Returns:
      Dictionary containing FASTA result, from polymer entity id to the
      `FastaSequence` object associated with that entity.
    """

    if verbosity:
        print("Querying RCSB for the '{}' FASTA file.".format(rcsb_id))
    response = requests.get(FASTA_BASE_URL + rcsb_id)

    if not response.ok:
        warnings.warn("It appears request failed with:" + response.text)
        response.raise_for_status()

    return _parse_fasta_text_to_list(response.text)
