#!/usr/bin/env python
import os
import sys
import argparse
import logging
import rpy2.robjects as robjects

import utils

rscript = ''
R = robjects.r


def run_rscript(command=None):
    """Run R command, log it, append to rscript"""
    global rscript
    if not command:
        return
    logging.debug(command)
    rscript += '{}\n'.format(command)
    output = R(command)
    return output


def prep_riboseqr_input(sam_file, output_file):
    """Generate input file for riboSeqR from SAM format file."""
    with open(output_file, 'w') as f:
        for line in open(sam_file):
            if line.startswith('@'):
                continue
            line = line.split()
            flag = line[1]

            if flag != '0':
                continue
            # make start 0-indexed, sam alignments are 1-indexed
            start = int(line[3]) - 1
            (name, sequence) = (line[2], line[9])
            f.write('"+"\t"{0}"\t{1}\t"{2}"\n'.format(name, start, sequence))


def batch_process(sam_files, seq_type, output_path):
    """Batch process the conversion of SAM format files -> riboSeqR format
    input files.

    Files are saved with file names corresponding to their sequence type.

    """
    outputs = []
    prefix = '{}'
    if seq_type == 'riboseq':
        prefix = 'RiboSeq file {}'
    elif seq_type == 'rnaseq':
        prefix = 'RNASeq file {}'

    for count, fname in enumerate(sam_files):
        count += 1
        out_file = os.path.join(output_path, prefix.format(count))
        logging.debug('Processing: {}'.format(fname))
        logging.debug('Writing output to: {}'.format(out_file))
        prep_riboseqr_input(fname, out_file)
        outputs.append(out_file)
    return outputs


def generate_ribodata(ribo_files='', rna_files='', replicate_names='',
                      seqnames='', rdata_save='Prepare.rda', sam_format=True,
                      html_file='Prepare-report.html', output_path=os.getcwd()):
    """Prepares Ribo and RNA seq data in the format required for riboSeqR. Calls
    the readRibodata function of riboSeqR and saves the result objects in an
    R data file which can be used as input for the next step.

    """
    input_ribo_files = utils.process_args(ribo_files, ret_mode='list')
    logging.debug('Found {} Ribo-Seq files'.format(len(input_ribo_files)))
    logging.debug(input_ribo_files)

    input_rna_files = []
    if rna_files:
        input_rna_files = utils.process_args(rna_files, ret_mode='list')
        logging.debug('Found {} RNA-Seq files'.format(len(input_rna_files)))
        logging.debug(input_rna_files)

    replicates = utils.process_args(replicate_names, ret_mode='charvector')
    logging.debug('Replicates: {}\n'.format(replicates))

    if sam_format:
        ribo_seq_files = batch_process(input_ribo_files, 'riboseq', output_path)
    else:
        ribo_seq_files = input_ribo_files

    html = '<h2>Prepare riboSeqR input - results</h2><hr>'
    if len(ribo_seq_files):
        html += '<h4>Generated riboSeqR format input files ' \
                '<em>(RiboSeq)</em></h4><p>'
        for fname in ribo_seq_files:
            html += '<a href="{0}">{0}</a><br>'.format(
                os.path.basename(fname))
        html += '</p>'

    rna_seq_files = []
    if len(input_rna_files):
        if sam_format:
            rna_seq_files = batch_process(
                input_rna_files, 'rnaseq', output_path)
        else:
            rna_seq_files = input_rna_files

    if len(rna_seq_files):
        html += ('<h4>Generated riboSeqR format input files '
                 '<em>(RNASeq)</em></h4><p>')
        for fname in rna_seq_files:
            html += '<a href="{0}">{0}</a><br>'.format(
                os.path.basename(fname))
        html += '</p>'

    input_seqnames = utils.process_args(seqnames, ret_mode='charvector')
    options = {'ribo_seq_files': 'c({})'.format(str(ribo_seq_files)[1:-1]),
               'rna_seq_files': 'c({})'.format(str(rna_seq_files)[1:-1]),
               'input_replicates': replicates,
               'input_seqnames': input_seqnames}

    logging.debug('{}'.format(R('sessionInfo()')))

    script = ''
    cmd = 'suppressMessages(library(riboSeqR))'
    run_rscript(cmd)
    script += '{}\n'.format(cmd)

    if len(rna_seq_files):
        cmd_args = ('riboFiles={ribo_seq_files}, '
                    'rnaFiles={rna_seq_files}'.format(**options))
    else:
        cmd_args = 'riboFiles={ribo_seq_files}'.format(**options)

    if input_seqnames:
        cmd_args += ', seqnames={input_seqnames}'.format(**options)
    if replicates:
        cmd_args += ', replicates={input_replicates}'.format(**options)
    else:
        cmd_args += ', replicates=c("")'
    cmd = 'riboDat <- readRibodata({0})'.format(cmd_args)
    run_rscript(cmd)
    script += '{}\n'.format(cmd)

    ribo_data = R['riboDat']
    logging.debug('riboDat \n{}\n'.format(ribo_data))
    cmd = 'save("riboDat", file="{}", compress=FALSE)'.format(rdata_save)
    run_rscript(cmd)
    script += '{}\n'.format(cmd)

    msg = '\n{:#^80}\n{}\n{:#^80}\n'.format(
        ' R script for this session ', script, ' End R script ')
    logging.debug(msg)

    with open(os.path.join(output_path, 'prepare.R'), 'w') as r:
        r.write(script)

    html += ('<h4>R script for this session</h4>'
             '<p><a href="prepare.R">prepare.R</a></p>'
             '<p>Next step: <em>Triplet periodicity</em></p>')

    with open(html_file, 'w') as f:
        f.write(html)

    return ribo_data


if __name__ == '__main__':

    description = (
        'Prepare riboSeqR input file from SAM format RNA/Ribo-Seq alignment.')
    parser = argparse.ArgumentParser(description=description)

    # required arguments
    flags = parser.add_argument_group('required arguments')
    flags.add_argument('--ribo_files', required=True,
                       help='List of Ribo-Seq files, comma-separated')
    # optional arguments
    parser.add_argument('--rna_files',
                        help='List of RNA-Seq files. Comma-separated')
    parser.add_argument('--replicate_names',
                        help='Replicate names, comma-separated')
    parser.add_argument('--seqnames',
                        help='Transcript (seqname) names to be read')
    parser.add_argument('--rdata_save',
                        help='File to write RData to (default: %(default)s)',
                        default='Prepare.rda')
    parser.add_argument(
        '--sam_format',
        help='Flag. Input is in SAM format', action='store_true')
    parser.add_argument('--html_file', help='Output file for results (HTML)')
    parser.add_argument('--output_path',
                        help='Files are saved in this directory')
    parser.add_argument('--debug', help='Flag. Produce debug output',
                        action='store_true')
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(format='%(module)s: %(levelname)s - %(message)s',
                            level=logging.DEBUG, stream=sys.stdout)
        logging.debug('Supplied Arguments: {}'.format(vars(args)))

    if not os.path.exists(args.output_path):
        os.mkdir(args.output_path)

    generate_ribodata(
        ribo_files=args.ribo_files, rna_files=args.rna_files,
        replicate_names=args.replicate_names, seqnames=args.seqnames,
        rdata_save=args.rdata_save,
        sam_format=args.sam_format, html_file=args.html_file,
        output_path=args.output_path
    )
    logging.debug('Done')
