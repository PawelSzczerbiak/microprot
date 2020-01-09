import os
import click
from microprot.scripts import process_fasta
from skbio import io


@click.command()
@click.option('--inp_fp', '-f', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='File path to a multi-sequence FASTA file.')
@click.option('--inp_identifier', '-d', required=True,
              help='Input identifier containing separate sequence ids and ranges '
                   'of sequence ids separated by comma e.g. 1-10,20,22-30 etc.')
@click.option('--microprot_inp', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Input directory where individual files from inp_fp '
                   'will be placed.')
@click.option('--microprot_out', '-o', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output directory path where processed_sequences.fasta '
                   'file will be created.')
def parse_inputs(inp_fp, inp_identifier, microprot_inp, microprot_out):
    """ Parse multi-sequence FASTA file into single-sequence, remove any
    problematic characters from the name and add information to
    `processed_sequences.fasta` file and `mapping.tsv` file.

    Returns
    -------
    SEQ_ids : list of str
        list of sequence ids picked from the inp_fp
    """
    for _dir in [microprot_inp, microprot_out]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    inp_ids = set()  # remove duplicates
    ranges = inp_identifier.split(',')
    for ran in ranges:
        el = ran.split('-')
        if len(el) == 1:  # separate id
            inp_ids.add(int(el[0]))
        else:  # proper range
            [inp_ids.add(i) for i in range(int(el[0]), int(el[1]) + 1)]
    inp_ids = sorted(inp_ids)  # ensure right order

    SEQS = process_fasta.extract_sequences(inp_fp, identifiers=inp_ids)
    SEQ_ids = []
    processed_fh = open('%s/%s' % (microprot_out,
                                   'processed_sequences.fasta'), 'a')
    mapping_fh = open('%s/%s' % (microprot_out,
                                 'mapping.tsv'), 'a')
    for i, SEQ in enumerate(SEQS):
        _seq = SEQ.metadata['id']
        _seq = _seq.replace('/', '_')
        _seq = _seq.replace('\\', '_')
        _seq = _seq.replace('|', '_')
        SEQ_ids.append(_seq)
        SEQ.metadata['id'] = _seq
        io.write(SEQ, format='fasta', into='%s/%s.fasta' % (microprot_inp,
                                                            _seq))
        io.write(SEQ, format='fasta',
                 into=processed_fh)
        mapping_fh.write("%i\t%s\n" % (inp_ids[i], _seq))
    processed_fh.close()
    mapping_fh.close()
    return SEQ_ids


if __name__ == "__main__":
    parse_inputs()
