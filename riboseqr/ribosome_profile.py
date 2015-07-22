#!/usr/bin/env python
import os
import sys
import glob
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
    msg = R(command)


def plot_transcript(rdata_load='Metagene.rda', transcript_name='',
                    transcript_length='27', transcript_cap='',
                    html_file='Plot-ribosome-profile.html',
                    output_path=os.getcwd()):
    """Plot ribosome profile for a given transcript. """
    options = {}
    for key, value, rtype, rmode in (
            ('transcript_name', transcript_name, 'str', None),
            ('transcript_length', transcript_length, 'int', 'charvector'),
            ('transcript_cap', transcript_cap, 'int', None)):
        options[key] = utils.process_args(value, ret_type=rtype, ret_mode=rmode)

    run_rscript('suppressMessages(library(riboSeqR))')
    run_rscript('load("{}")'.format(rdata_load))

    html = """<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 3.2//EN">
    <html>
    <head>
    <title>Ribosome Profile Plot - Report</title>
    </head>
    <body>
    """
    html += '<h2>Plot ribosome profile - results</h2>\n<hr>\n'
    if len(transcript_name):
        cmd_args = (
            '"{transcript_name}", main="{transcript_name}",'
            'coordinates=ffCs@CDS, riboData=riboDat,'
            'length={transcript_length}'.format(**options))
        if transcript_cap:
            cmd_args += ', cap={transcript_cap}'.format(**options)
        plot_file = os.path.join(output_path, 'Ribosome-profile-plot')
        for fmat in ('pdf', 'png'):
            if fmat == 'png':
                cmd = 'png(file="{}_%1d.png", type="cairo")'.format(plot_file)
            else:
                cmd = 'pdf(file="{}.pdf")'.format(plot_file)
            run_rscript(cmd)
            cmd = 'plotTranscript({})'.format(cmd_args)
            run_rscript(cmd)
            run_rscript('dev.off()')

        html += ('<p>Selected ribosome footprint length: '
                 '<strong>{0}</strong>\n'.format(transcript_length))

        for image in sorted(glob.glob('{}_*.png'.format(plot_file))):
            html += '<p><img border="1" src="{0}" alt="{0}"></p>\n'.format(
                os.path.basename(image))
        html += '<p><a href="Ribosome-profile-plot.pdf">PDF version</a></p>\n'
    else:
        msg = 'No transcript name was provided. Did not generate plot.'
        html += '<p>{}</p>'.format(msg)
        logging.debug(msg)

    logging.debug('\n{:#^80}\n{}\n{:#^80}\n'.format(
        ' R script for this session ', rscript, ' End R script '))

    with open(os.path.join(output_path, 'ribosome-profile.R'), 'w') as r:
        r.write(rscript)

    html += ('<h4>R script for this session</h4>\n'
             '<p><a href="ribosome-profile.R">ribosome-profile.R</a></p>\n'
             '</body>\n</html>\n')

    with open(html_file, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot Ribosome profile')

    # required arguments
    flags = parser.add_argument_group('required arguments')
    flags.add_argument('--rdata_load', required=True,
                       help='Saved riboSeqR data from Step 2')
    flags.add_argument('--transcript_name', required=True,
                       help='Name of the transcript to be plotted')
    flags.add_argument(
        '--transcript_length', required=True,
        help='Size class of ribosome footprint data to be plotted',
        default='27')
    flags.add_argument(
        '--transcript_cap', required=True,
        help=('Cap on the largest value that will be plotted as an abundance '
              'of the ribosome footprint data'))
    parser.add_argument('--html_file', help='HTML file with reports')
    parser.add_argument('--output_path', help='Directory to save output files')
    parser.add_argument('--debug', help='Produce debug output',
                        action='store_true')

    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(format='%(levelname)s - %(message)s',
                            level=logging.DEBUG, stream=sys.stdout)
        logging.debug('Supplied Arguments\n{}\n'.format(vars(args)))

    if not os.path.exists(args.output_path):
        os.mkdir(args.output_path)

    plot_transcript(rdata_load=args.rdata_load,
                    transcript_name=args.transcript_name,
                    transcript_length=args.transcript_length,
                    transcript_cap=args.transcript_cap,
                    html_file=args.html_file, output_path=args.output_path)
    logging.debug('Done!')
