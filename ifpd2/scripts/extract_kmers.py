"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from Bio import SeqIO  # type: ignore
from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
import click  # type: ignore
from ifpd2.const import CONTEXT_SETTINGS
import logging
from os.path import join as path_join
from os.path import basename, normpath, splitext


@click.command(
    name="extract_kmers",
    context_settings=CONTEXT_SETTINGS,
    help="Generate oligonucleotides K-mers from FASTA",
)
@click.argument("input_path", metavar="INPUT_FASTA", type=click.Path(exists=True))
@click.argument("output_path", metavar="OUTPUT_DIRECTORY", type=click.Path(exists=True))
@click.argument("kmer_size", metavar="KMER_LENGTH", type=click.INT)
def main(input_path: str, output_path: str, kmer_size: int) -> None:
    seqRec = SeqIO.parse(input_path, "fasta")

    oligos_list = []
    for record in seqRec:
        record = record
        for i in range(len(record) - kmer_size + 1):
            oligos_list.append(
                SeqRecord(
                    Seq(str(record.seq)[slice(i, i + kmer_size)]).reverse_complement(),
                    id=f"{record.id}|{i+1}:{i+kmer_size+1}",
                    name="",
                    description="",
                )
            )
    logging.info(f"Extracted {len(oligos_list)} sequences")

    valid_oligos = [
        oligo
        for oligo in oligos_list
        if (oligo.seq.count("C") + oligo.seq.count("G")) / kmer_size >= 0.35
        and (oligo.seq.count("C") + oligo.seq.count("G")) / kmer_size <= 0.85
        and "N" not in oligo.seq
    ]

    print("Sequences with correct GC content:", len(valid_oligos))

    base, _ = splitext(basename(input_path))
    SeqIO.write(
        valid_oligos,
        path_join(normpath(output_path), f"{base}.GC35to85_RevCompl.fa"),
        "fasta",
    )

    oligos_rc = [
        oligo.reverse_complement(id=True, name="", description="")
        for oligo in valid_oligos
    ]
    SeqIO.write(
        oligos_rc,
        path_join(normpath(output_path), f"{base}.GC35to85_Reference.fa"),
        "fasta",
    )
