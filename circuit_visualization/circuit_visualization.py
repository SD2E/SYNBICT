import logging
import argparse
import os
import sys

from sbol import *
import math
import dnaplotlib as dpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Polygon, Ellipse, Wedge, Circle, PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.patheffects import Stroke
import matplotlib.patches as patches

from sequences_to_features import FeatureLibrary
from features_to_circuits import CircuitLibrary

def curved_activation(ax, type, num, from_part, to_part, scale, linewidth, arc_height_index, opts):
    color = (0.0,0.0,0.0)
    arcHeightStart = 10
    arcHeightEnd = 10

    if opts != None:
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'arc_height_start' in list(opts.keys()):
            arcHeightStart = opts['arc_height_start']
        if 'arc_height_end' in list(opts.keys()):
            arcHeightEnd = opts['arc_height_end']
        if 'rad' in list(opts.keys()):
            rad = opts['rad']

    start = (from_part['start'] + from_part['end'])/2
    end = (to_part['start'] + to_part['end'])/2

    if start > end:
        arcHeightStart = -arcHeightStart
        arcHeightEnd = -arcHeightEnd

    ax.annotate('', (end, arcHeightEnd), (start, arcHeightStart), ha="right", va="center", size=8,
        arrowprops=dict(arrowstyle='->', connectionstyle="arc3,rad=" + str(rad), lw=linewidth, color=color))

def curved_repression(ax, type, num, from_part, to_part, scale, linewidth, arc_height_index, opts):
    color = (0.0,0.0,0.0)
    arcHeightStart = 10
    arcHeightEnd = 10
    rad = -0.4

    if opts != None:
        if 'linewidth' in list(opts.keys()):
            linewidth = opts['linewidth']
        if 'color' in list(opts.keys()):
            color = opts['color']
        if 'arc_height_start' in list(opts.keys()):
            arcHeightStart = opts['arc_height_start']
        if 'arc_height_end' in list(opts.keys()):
            arcHeightEnd = opts['arc_height_end']
        if 'rad' in list(opts.keys()):
            rad = opts['rad']

    start = (from_part['start'] + from_part['end'])/2
    end = (to_part['start'] + to_part['end'])/2

    if start > end:
        arcHeightStart = -arcHeightStart
        arcHeightEnd = -arcHeightEnd

    ax.annotate('', (end, arcHeightEnd), (start, arcHeightStart), ha="right", va="center", size=8,
        arrowprops=dict(arrowstyle='-', connectionstyle="arc3,rad=" + str(rad), lw=linewidth, color=color))

    arrowhead_length = 4

    delta = arrowhead_length/math.sqrt(1 + 4*math.pow(rad, 2))

    if rad < 0:
        delta_y = delta
        delta_x = -2*rad*delta
    else:
        delta_y = -delta
        delta_x = 2*rad*delta

    head_x = end - 1.1*arrowhead_length
    head_y = arcHeightEnd + 0.95*arrowhead_length

    line_rep = Line2D([head_x - delta_x/2, head_x + delta_x/2], [head_y - delta_y/2, head_y + delta_y/2], 
        linewidth=linewidth, color=color, zorder=12, linestyle='-')

    ax.add_line(line_rep)

def load_sbol(sbol_file):
    print('Loading ' + sbol_file)

    doc = Document()
    doc.read(sbol_file)

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')

    return doc

class CircuitVisualizer():

    SO_DICT = {
        SO_TERMINATOR: 'Terminator',
        SO_CDS: 'CDS',
        SO_PROMOTER: 'Promoter',
        SO_RBS: 'RBS'
    }

    COLOR = {
        'red': (0.95, 0.30, 0.25),
        'green': (0.38, 0.82, 0.32),
        'blue': (0.38, 0.65, 0.87),
        'orange': (1.00, 0.75, 0.17),
        'purple': (0.55, 0.35, 0.64)
    }

    def __init__(self):
        pass

    @classmethod
    def visualize(cls, circuit_library, feature_library, min_features):
        lw = 0.75

        dr = dpl.DNARenderer()

        part_renderers = dr.SBOL_part_renderers()

        reg_renderers = dr.std_reg_renderers()
        reg_renderers['CurvedRepression'] = curved_repression
        reg_renderers['CurvedActivation'] = curved_activation

        k = 0

        for circuit in circuit_library.circuits:
            if len(circuit.features) > min_features:
                print('Visualizing ' + circuit.identity)

                k = k + 1

                designs = []

                part_dict = {}

                for feature in circuit.features:
                    comp_definition = feature_library.get_definition(feature.identity)

                    if comp_definition is not None and BIOPAX_DNA in comp_definition.types:
                        design = []

                        located = []

                        for seq_anno in comp_definition.sequenceAnnotations:
                            if seq_anno.component is not None and len(seq_anno.locations) == 1 and seq_anno.locations[0].getTypeURI() == SBOL_RANGE:
                                sub_comp = comp_definition.components.get(seq_anno.component)

                                range_location = seq_anno.locations.getRange()

                                located.append((range_location.start, range_location.end, range_location.orientation, sub_comp.definition,
                                    sub_comp.identity))

                        located.sort()

                        max_span = 1

                        for i in range(1, len(located)):
                            span = located[i][0] - located[i - 1][1]

                            if span > max_span:
                                max_span = span

                        for i in range(0, len(located)):
                            if i > 0:
                                span = located[i][0] - located[i - 1][1]

                                if span > 0:
                                    normal_span = int(span/max_span*60)
                                else:
                                    normal_span = 0

                                if normal_span > 0:
                                    design.append({'type': 'EmptySpace', 'fwd': True, 'opts': {'x_extent': normal_span}})

                            feature_definition = feature_library.get_definition(located[i][3])

                            part = {'opts': {}}

                            for role in feature_definition.roles:
                                if role in cls.SO_DICT:
                                    part['type'] = cls.SO_DICT[role]
                                elif 'type' not in part:
                                    part['type'] = 'UserDefined'

                            if part['type'] == 'CDS':
                                part['opts']['arrowhead_height'] = 0

                            part['name'] = located[i][4]

                            if part['type'] == 'CDS' or part['type'] == 'promoter':
                                if feature_definition.name is None:
                                    part['opts']['label'] = feature_definition.displayId
                                else:
                                    part['opts']['label'] = feature_definition.name

                                part['opts']['label_y_offset'] = -30
                                part['opts']['label_size'] = 4
                                if part['type'] == 'CDS':
                                    part['opts']['label_style'] = 'italic'

                            part['fwd'] = (located[i][2] == SBOL_ORIENTATION_INLINE)

                            if located[i][3] not in part_dict:
                                part_dict[located[i][3]] = []

                            part_dict[located[i][3]].append(part)

                            design.append(part)

                        designs.append(design)

                arcs = []
                
                for template in circuit.get_feature_identities():
                    for activated in circuit_library.get_activated_by_template(template):
                        for source in part_dict[template]:
                            for sink in part_dict[activated]:
                                arcs.append({'type':'CurvedActivation', 'from_part': source, 'to_part': sink, 'opts': {'color': cls.COLOR['green'],
                                    'linewidth':lw, 'arc_height_start': 15, 'arc_height_end': 15, 'rad': -0.4}})

                    for repressed in circuit_library.get_repressed_by_template(template):
                        for source in part_dict[template]:
                            for sink in part_dict[repressed]:
                                arcs.append({'type':'CurvedRepression', 'from_part': source, 'to_part': sink, 'opts': {'color': cls.COLOR['red'],
                                    'linewidth':lw, 'arc_height_start': 15, 'arc_height_end': 15, 'rad': -0.4}})

                fig = plt.figure(figsize=(2.2,1.8))
                gs = gridspec.GridSpec(len(designs), 1)

                dna_axes = []

                for i in range(0, len(designs)):
                    dna_axes.append(plt.subplot(gs[i]))

                for i in range(0, len(designs)):
                    dna_axis = plt.subplot(gs[i])
                    start, end = dr.renderDNA(dna_axis, designs[i], part_renderers, regs=arcs, reg_renderers=reg_renderers)
                    # start, end = dr.renderDNA(ax_dna1, design1, part_renderers, regs=reg1, reg_renderers=reg_renderers)
                    dna_axis.set_xlim([start, end])
                    dna_axis.set_ylim([-25,25])
                    dna_axis.set_aspect('equal')
                    dna_axis.set_xticks([])
                    dna_axis.set_yticks([])
                    dna_axis.axis('off')

                plt.subplots_adjust(hspace=0.01, left=0.05, right=0.95, top=0.92, bottom=0.01)

                fig.savefig('test_circuit' + str(k) + '.pdf', transparent=True)
                fig.savefig('test_circuit' + str(k) + '.png', dpi=300)

                plt.close('all')

                logging.info('Finished visualizing %s.\n', circuit.identity)
            else:
                logging.warning('Failed to visualize %s.\n', circuit.identity)

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--circuit_file')
    parser.add_argument('-f', '--feature_files', nargs='*', default=[])
    parser.add_argument('-l', '--curation_log', nargs='?', default='')
    parser.add_argument('-m', '--min_features', nargs='?', default=2)
    parser.add_argument('-v', '--validate', action='store_true')
    
    args = parser.parse_args(args)

    if len(args.curation_log) > 0:
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(level=logging.INFO, filename=args.curation_log, filemode='w',
                            format='%(levelname)s : %(message)s')

    circuit_doc = load_sbol(args.circuit_file)

    circuit_library = CircuitLibrary([circuit_doc])

    feature_docs = []

    for feature_file in args.feature_files:
        feature_docs.append(load_sbol(feature_file))
    feature_library = FeatureLibrary(feature_docs)

    CircuitVisualizer.visualize(circuit_library, feature_library, args.min_features)

    print('Finished visualizing.')

if __name__ == '__main__':
    main()
