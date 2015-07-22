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


def find_periodicity(
        rdata_load='Prepare.rda', start_codons='ATG', stop_codons='TAG,TAA,TGA',
        fasta_file=None, include_lengths='25:30', analyze_plot_lengths='26:30',
        text_legend='Frame 0, Frame 1, Frame 2', rdata_save='Periodicity.rda',
        html_file='Periodicity-report.html', output_path=os.getcwd()):
    """Plot triplet periodicity from prepared R data file. """
    logging.debug('{}'.format(R('sessionInfo()')))
    cmd = 'suppressMessages(library(riboSeqR))'
    run_rscript(cmd)

    logging.debug('Loading saved R data file')
    cmd = 'load("{}")'.format(rdata_load)
    run_rscript(cmd)

    # R("""options(showTailLines=Inf)""")
    starts, stops = (utils.process_args(start_codons, ret_mode='charvector'),
                     utils.process_args(stop_codons, ret_mode='charvector'))

    cmd = ('fastaCDS <- findCDS(fastaFile={0!r}, startCodon={1}, '
           'stopCodon={2})'.format(fasta_file, starts, stops))
    run_rscript(cmd)

    logging.debug('Potential coding sequences using start codon (ATG) and '
                  'stop codons TAG, TAA, TGA')
    logging.debug('{}\n'.format(R['fastaCDS']))

    cmd = """fCs <- frameCounting(riboDat, fastaCDS, lengths={0})
    fS <- readingFrame(rC=fCs, lengths={1}); fS""".\
        format(include_lengths, analyze_plot_lengths)
    run_rscript(cmd)

    logging.debug('riboDat \n{}\n'.format(R['riboDat']))
    logging.debug('fCs\n{0}\n'.format(R['fCs']))
    logging.debug('Reading frames for each n-mer\n{}'.format(R['fS']))

    legend = utils.process_args(text_legend, ret_mode='charvector')

    for fmat in ('pdf', 'png'):
        if fmat == 'png':
            cmd = '{0}(file="{1}", type="cairo")'
        else:
            cmd = '{0}(file="{1}")'
        run_rscript(cmd.format(fmat, os.path.join(
            output_path, '{0}.{1}'.format('Periodicity-plot', fmat))))
        run_rscript('plotFS(fS, legend.text = {0})'.format(legend))
        run_rscript('dev.off()')

    run_rscript('save("fCs", "fS", "riboDat", "fastaCDS", '
                'file="{}", compress=FALSE)'.format(rdata_save))

    html = '<h2>Triplet periodicity - results</h2><hr>'
    html += ('<h4>Results of reading frame analysis</h4>'
             '<pre>{}</pre><br>'.format(R['fS']))
    html += ('<p>Lengths used for reading frame analysis - <code>{0}</code>'
             '<br>Lengths selected for the plot - <code>{1}</code>'
             '</p>'.format(include_lengths, analyze_plot_lengths))
    html += ('<p><img src="Periodicity-plot.png" border="1" '
             'alt="Triplet periodicity plot" />'
             '<br><a href="Periodicity-plot.pdf">PDF version</a></p>')

    logging.debug('\n{:#^80}\n{}\n{:#^80}\n'.format(
        ' R script for this session ', rscript, ' End R script '))

    with open(os.path.join(output_path, 'periodicity.R'), 'w') as r:
        r.write(rscript)

    html += ('<h4>R script for this session</h4>'
             '<p><a href="periodicity.R">periodicity.R</a></p>'
             '<p>Next step: <em>Metagene analysis</em></p>')

    with open(html_file, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description='Plot triplet periodicity for different read lengths.')
    
    # required arguments
    flags = parser.add_argument_group('required arguments')
    flags.add_argument(
        '--rdata_load', required=True,
        help='riboSeqR data from input preparation step (Prepare.rda)')
    flags.add_argument(
        '--fasta_file', required=True,
        help='FASTA file of the reference transcriptome')
    
    # optional arguments
    parser.add_argument(
        '--start_codons',
        help='Start codon(s) to use. (default: %(default)s)', default='ATG')
    parser.add_argument(
        '--stop_codons', help='Stop codon(s) to use. (default: %(default)s)',
        default='TAG, TAA, TGA')
    parser.add_argument(
        '--include_lengths',
        help='Lengths of ribosome footprints to be included '
             '(default: %(default)s)', default='25:30')
    parser.add_argument(
        '--analyze_plot_lengths',
        help='Lenghts of reads to be analyzed for frame-shift and plotting '
             '(default: %(default)s)', default='26:30')
    parser.add_argument(
        '--text_legend',
        help='Text for legend used in the plot (default: %(default)s)',
        default='Frame 0, Frame 1, Frame 2')
    parser.add_argument(
        '--rdata_save', help='File to write RData to (default: %(default)s)',
        default='Periodicity.rda')
    parser.add_argument('--html_file', help='Output file for results (HTML)')
    parser.add_argument('--output_path',
                        help='Files are saved in this directory')
    parser.add_argument(
        '--debug', help='Produce debug output', action='store_true')
    
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(format='%(module)s: %(levelname)s - %(message)s',
                            level=logging.DEBUG, stream=sys.stdout)
        logging.debug('Supplied Arguments\n{}\n'.format(vars(args)))

    if not os.path.exists(args.output_path):
        os.mkdir(args.output_path)

    find_periodicity(
        rdata_load=args.rdata_load, start_codons=args.start_codons,
        stop_codons=args.stop_codons, fasta_file=args.fasta_file,
        include_lengths=args.include_lengths,
        analyze_plot_lengths=args.analyze_plot_lengths,
        text_legend=args.text_legend,
        rdata_save=args.rdata_save, html_file=args.html_file,
        output_path=args.output_path)
logging.debug("Done!")
