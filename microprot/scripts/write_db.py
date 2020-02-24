import os
import textwrap
from datetime import datetime
import click

from microprot.scripts import process_fasta


def msa_size(msa_fp):
    """ Determine size of an MSA
    Parameters
    ----------
    msa_fp : str
        File path to an MSA file (a3m or root of file name)

    Returns
    -------
    msa_size : int
        size of an MSA
    """
    msa_dir, msa_ext = os.path.splitext(os.path.abspath(msa_fp))
    if msa_ext != '.a3m':
        msa_ext = '.a3m'
        msa_fp = ''.join([msa_dir, msa_ext])

    with open(msa_fp, 'r') as f:
        lines = f.readlines()
    msa_size = sum(1 for line in lines if line.startswith('>')) - 1
    return msa_size


@click.command()
@click.option('--fname', '-f', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='File path to a multi-sequence FASTA file.')
@click.option('--step', '-s', required=True,
              help='Processing step information (e.g. CM, Pfam).')
@click.option('--version', '-v', required=False, default="1",
              help='Processing version.')
@click.option('--db_fp', '-o', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output database information. Sequence with header '
                   ' will be appended  to `db_fp` and header with '
                   'processing information do `db_fp.index`.')
def write_db(fname, step, version, db_fp):
    """Append protein name and other information into the
    sequence DB Parameters"""

    prots = process_fasta.extract_sequences(fname)
    for prot in prots:
        prot_name = prot.metadata['id']
        timestamp = str(datetime.now()).split('.')[0]
        msa_size_len = msa_size(fname)

        # > protein_name # source # msa_size # commit_no # timestamp
        append_idx = '>%s # %s # %i # %s # %s\n' % (prot_name,
                                                    step,
                                                    msa_size_len,
                                                    version,
                                                    timestamp
                                                    )
        with open('%s.index' % db_fp, 'a') as f:
            f.write(append_idx)

        append_seq = '>%s\n%s\n' % (prot_name, textwrap.fill(str(prot[:]), 70))
        with open(db_fp, 'a') as f:
            f.write(append_seq)


if __name__ == "__main__":
    write_db()
