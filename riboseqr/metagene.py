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
    output = R(command)
    return output


def do_analysis(
        rdata_load='Periodicity.rda', selected_lengths='27',
        selected_frames='', hit_mean='10', unique_hit_mean='1',
        ratio_check='TRUE', min5p='-20', max5p='200', min3p='-200', max3p='20',
        cap='', plot_title='', plot_lengths='27', rdata_save='Metagene.rda',
        html_file='Metagene-report.html', output_path=os.getcwd()):
    """Metagene analysis from saved periodicity R data file. """
    run_rscript('suppressMessages(library(riboSeqR))')
    run_rscript('load("{}")'.format(rdata_load))

    logging.debug('fS\n{}\nfCs\n{}\n'.format(R['fS'], R['fCs']))
    options = {}
    for key, value, rtype, rmode in (
            ('lengths', selected_lengths, 'int', 'charvector'),
            ('frames', selected_frames, 'int', 'listvector'),
            ('hit_mean', hit_mean, 'int', None),
            ('unique_hit_mean', unique_hit_mean, 'int', None),
            ('ratio_check', ratio_check, 'bool', None),
            ('min5p', min5p, 'int', None), ('max5p', max5p, 'int', None),
            ('min3p', min3p, 'int', None), ('max3p', max3p, 'int', None),
            ('cap', cap, 'int', None),
            ('plot_title', plot_title, 'str', 'charvector'),
            ('plot_lengths', plot_lengths, 'int', 'list')):
            options[key] = utils.process_args(
                value, ret_type=rtype, ret_mode=rmode)

    cmd_args = """fCs, lengths={lengths},
    frames={frames}, hitMean={hit_mean},
    unqhitMean={unique_hit_mean}, fS=fS""".format(**options)

    if ratio_check == 'TRUE':
        cmd_args += ', ratioCheck = TRUE'

    run_rscript('ffCs <- filterHits({})'.format(cmd_args))
    logging.debug("ffCs\n{}\n".format(R['ffCs']))

    cds_args = ('coordinates=ffCs@CDS, riboDat=riboDat, min5p={min5p}, '
                'max5p={max5p}, min3p={min3p}, max3p={max3p}'.format(**options))

    if options['cap']:
        cds_args += ', cap={cap}'.format(**options)

    if options['plot_title']:
        cds_args += ', main={plot_title}'.format(**options)

    html = """<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 3.2//EN">
    <html>
    <head>
    <title>Metagene Analysis - Report</title>
    </head>
    <body>
    """
    html += '<h2>Metagene analysis - results</h2>\n<hr>\n'
    html += ('<p>\nLengths of footprints used in analysis - <strong>'
             '<code>{0}</code></strong><br>\nLengths of footprints '
             'selected for the plot - <strong><code>{1}</code></strong>'
             '\n</p>\n'.format(selected_lengths, plot_lengths))
    for count, length in enumerate(options['plot_lengths']):
        count += 1
        html += '<h3>Length: {0}</h3>\n'.format(length)
        plot_file = os.path.join(output_path,
                                 'Metagene-analysis-plot{0}'.format(count))
        for fmat in ('pdf', 'png'):
            if fmat == 'png':
                cmd = 'png(file="{0}_%1d.png", type="cairo")'
            else:
                cmd = 'pdf(file="{0}.pdf")'
            run_rscript(cmd.format(plot_file))
            run_rscript('plotCDS({0},{1})'.format(
                cds_args, 'lengths={}'.format(length)))
            run_rscript('dev.off()')
        for image in sorted(
                glob.glob('{}*.png'.format(plot_file))):
            html += '<p><img border="1" src="{0}" alt="{0}"></p>\n'.format(
                os.path.basename(image))
        html += '<p><a href="{0}.pdf">PDF version</a></p>\n'.format(
            os.path.basename(plot_file))
    run_rscript('save("ffCs", "riboDat", "fastaCDS", file="{}", '
                'compress=FALSE)'.format(rdata_save))

    logging.debug('\n{:#^80}\n{}\n{:#^80}\n'.format(
        ' R script for this session ', rscript, ' End R script '))

    with open(os.path.join(output_path, 'metagene.R'), 'w') as r:
        r.write(rscript)

    html += ('<h4>R script for this session</h4>\n'
             '<p><a href="metagene.R">metagene.R</a></p>\n'
             '<p>Next step: <em>Plot Ribosome profile</em></p>\n'
             '</body>\n</html>\n')

    with open(html_file, 'w') as f:
        f.write(html)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Metagene analysis')

    # required arguments
    flags = parser.add_argument_group('required arguments')
    flags.add_argument('--rdata_load', required=True,
                       help='Saved riboSeqR data from Periodicity step.')
    flags.add_argument('--selected_lengths', required=True,
                       help='Select frame lengths to filter. Comma-separated',
                       default='27')
    flags.add_argument(
        '--selected_frames', required=True,
        help='Select frames corresponding to frame lengths. Comma-separated')

    flags.add_argument(
        '--hit_mean', required=True,
        help='Mean number of hits within the replicate group for filtering',
        default='10')

    flags.add_argument(
        '--unique_hit_mean', required=True,
        help='Mean number of unique sequences within the replicate group '
             'for filtering', default='1')

    parser.add_argument(
        '--rdata_save', help='File to write R data to (default: %(default)s)',
        default='Metagene.rda')

    parser.add_argument(
        '--ratio_check',
        help='Check the ratios of the expected phase to maximal phase '
             'within the putative coding sequence (default: %(default)s)',
        choices=['TRUE', 'FALSE'], default='TRUE')

    parser.add_argument(
        '--plot_lengths',
        help='Length of footprints to be plotted. Multiple values should be '
             'comma-separated. In that case, multiple plots will be produced'
             '(default: %(default)s)', default='27')

    parser.add_argument(
        '--min5p',
        help='The distance upstream of the translation start to be plotted '
             '(default: %(default)s)', default='-20')

    parser.add_argument(
        '--max5p',
        help='The distance downstream of the translation start to be plotted '
             '(default: %(default)s)', default='200')

    parser.add_argument(
        '--min3p',
        help='The distance upstream of the translation end to be plotted '
             '(default: %(default)s)', default='-200')

    parser.add_argument(
        '--max3p',
        help='The distance downtream of the translation end to be plotted '
             '(default: %(default)s)', default='20')

    parser.add_argument(
        '--cap', help='If given, caps the height of plotted values '
                      '(default: %(default)s)')

    parser.add_argument('--plot_title', help='Title of the plot', default='')
    parser.add_argument('--html_file', help='HTML file with reports')
    parser.add_argument('--output_path', help='Directory to save output files')
    parser.add_argument(
        '--debug', help='Produce debug output', action='store_true')

    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(format='%(module)s: %(levelname)s - %(message)s',
                            level=logging.DEBUG, stream=sys.stdout)
        logging.debug('Supplied Arguments\n{}\n'.format(vars(args)))

    if not os.path.exists(args.output_path):
        os.mkdir(args.output_path)

    do_analysis(
        rdata_load=args.rdata_load, selected_lengths=args.selected_lengths,
        selected_frames=args.selected_frames, hit_mean=args.hit_mean,
        unique_hit_mean=args.unique_hit_mean, ratio_check=args.ratio_check,
        min5p=args.min5p, max5p=args.max5p, min3p=args.min3p, max3p=args.max3p,
        cap=args.cap, plot_title=args.plot_title,
        plot_lengths=args.plot_lengths, rdata_save=args.rdata_save,
        html_file=args.html_file, output_path=args.output_path)

    logging.debug('Done!')
