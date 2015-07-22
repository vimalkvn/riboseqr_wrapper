# -*- coding: utf-8 -*-
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


def get_counts(rdata_load='Metagene.rda', slice_lengths='27',
               frames='', group1=None, group2=None, num_counts=10,
               normalize='FALSE', html_file='Counts.html',
               output_path='counts'):
    options = {'slice_lengths': utils.process_args(
        slice_lengths, ret_type='int', ret_mode='charvector')}

    if not num_counts:
        num_counts = 10

    options['frames'] = utils.process_args(
        frames, ret_type='int', ret_mode='listvector')

    run_rscript('suppressMessages(library(riboSeqR))')
    run_rscript('load("{}")'.format(rdata_load))

    cmd_args = 'ffCs, lengths={slice_lengths}'.format(**options)
    if frames:
        cmd_args += ', frames={frames}'.format(**options)

    run_rscript("""riboCounts <- sliceCounts({})
    annotation <- as.data.frame(ffCs@CDS)
    rownames(riboCounts) <- annotation$seqnames
    colnames(riboCounts) <- names(riboDat@riboGR)""".format(cmd_args))

    run_rscript("""mrnaCounts <- rnaCounts(riboDat, ffCs@CDS)
    rownames(mrnaCounts) <- annotation$seqnames""")

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    html = '<h2>Differential Translation Analysis</h2><hr>'
    for count_name, file_name, legend in (
            ('riboCounts', 'RiboCounts.csv', 'Ribo-Seq counts'),
            ('mrnaCounts', 'RNACounts.csv', 'RNA-Seq counts'),
            ('tC', 'TopCounts.csv', 'baySeq topCounts')):
        if count_name == 'tC' and R['riboCounts'] and R['mrnaCounts']:
            run_rscript('suppressMessages(library(baySeq))')
            cmd = """pD <- new("countData", replicates=ffCs@replicates, \
            data=list(riboCounts, mrnaCounts), groups=list(NDT={0}, DT={1}), \
            annotation=as.data.frame(ffCs@CDS), \
            densityFunction=bbDensity)""".format(
                utils.process_args(
                    group1, ret_type='int', ret_mode='charvector'),
                utils.process_args(
                    group2, ret_type='str', ret_mode='charvector'))
            run_rscript(cmd)

            run_rscript('libsizes(pD) <- getLibsizes(pD)')
            run_rscript('pD <- getPriors(pD, cl=NULL)')
            run_rscript('pD <- getLikelihoods(pD, cl=NULL)')
            run_rscript('tC <- topCounts(pD, "DT", normaliseData={}, '
                        'number={})'.format(normalize, num_counts))

        if R[count_name]:
            html += '<h3>{}</h3>'.format(legend)
            output_file = os.path.join(output_path, file_name)
            run_rscript('write.csv({}, file="{}")'.format(
                count_name, output_file))
            with open(output_file) as g:
                html += '<table cellpadding="4" border="1">'
                header = g.readline()
                html += '<tr>'
                for item in header.strip().split(","):
                    html += '<th>{}</th>'.format(item.strip('"'))
                html += '</tr>'
                for line in g:
                    html += '<tr>'
                    data = line.strip().split(",")
                    # html += '<td>{}</td>'.format(data[0].strip('"'))
                    for item in data:
                        html += ('<td align="center"><code>{}</code>'
                                 '</td>'.format(item.strip('"')))
                    html += '</tr>'
                html += ('</table><p>Download: <a href="{0}">'
                         '{0}</a></p><hr>'.format(file_name))

    with open(os.path.join(output_path, 'counts.R'), 'w') as r:
        r.write(rscript)

    html += ('<h4>R script for this session</h4>'
             '<p>Download: <a href="counts.R">counts.R</a></p>')

    with open(html_file, 'w') as f:
        f.write(html)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get Ribo and RNA-Seq counts')

    # required arguments
    flags = parser.add_argument_group('required arguments')
    flags.add_argument('--rdata_load', required=True,
                       help='Saved riboSeqR data from Step 2')
    flags.add_argument('--slice_lengths',
                       help='Lengths of ribosome footprints to inform count '
                            'data', default='27,28', required=True)
    flags.add_argument('--group1', help='Group model for baySeq', required=True)
    flags.add_argument('--group2', help='Group model for baySeq', required=True)
    flags.add_argument('--frames',
                       help='Frames of ribosome footprints (relative to '
                            'coding start site). If omitted, all frames '
                            'are used.', required=True)

    # optional arguments
    parser.add_argument(
        '--num_counts', help='How many results to return? (topCounts)')
    parser.add_argument('--normalize', help='Normalize data?',
                        choices=['TRUE', 'FALSE'], default='FALSE')
    parser.add_argument('--html_file', help='HTML file with reports')
    parser.add_argument('--output_path', help='Directory to save output files')
    parser.add_argument(
        '--debug', help='Produce debug output', action='store_true')

    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(format='%(levelname)s - %(message)s',
                            level=logging.DEBUG, stream=sys.stdout)
        logging.debug('Supplied Arguments\n{}\n'.format(vars(args)))

    if not os.path.exists(args.output_path):
        os.mkdir(args.output_path)

    get_counts(rdata_load=args.rdata_load, slice_lengths=args.slice_lengths,
               frames=args.frames, group1=args.group1, group2=args.group2,
               num_counts=args.num_counts, normalize=args.normalize,
               html_file=args.html_file, output_path=args.output_path)
