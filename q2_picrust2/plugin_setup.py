from qiime2.plugin import (Plugin, Str, Choices, Int, Bool, Range, Float,
                           Citations)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Sequence
from q2_types.tree import Phylogeny, Rooted
import q2_picrust2

citations = Citations.load('citations.bib', package='q2_picrust2')

HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']

PLACEMENT_TOOLS = ['epa-ng', 'sepp']

plugin = Plugin(
    name='picrust2',
    version="2023.2",
    website='https://github.com/gavinmdouglas/q2-picrust2',
    package='q2_picrust2',
    description=('This QIIME 2 plugin wraps the default 16S PICRUSt2 pipeline to run '
                 'metagenome inference based on marker gene data. Currently '
                 'only unstratified output is supported.'),
    short_description='Predicts gene families and pathways from 16S sequences.',
    citations=[citations['Douglas2020NatureBiotech']]
)

plugin.methods.register_function(
    function=q2_picrust2.full_pipeline,

    inputs={'table': FeatureTable[Frequency],
            'seq': FeatureData[Sequence]},

    parameters={'threads': Int % Range(1, None),
                'hsp_method': Str % Choices(HSP_METHODS),
                'placement_tool': Str % Choices(PLACEMENT_TOOLS),
                'min_align': Float % Range(0.0, 1.0),
                'max_nsti': Float % Range(0.0, None),
                'edge_exponent': Float % Range(0.0, None),
                'skip_minpath': Bool,
                'no_gap_fill': Bool,
                'skip_norm': Bool,
                'highly_verbose': Bool},

    outputs=[('ko_metagenome', FeatureTable[Frequency]),
             ('ec_metagenome', FeatureTable[Frequency]),
             ('pathway_abundance', FeatureTable[Frequency])],

    input_descriptions={
        'table': ('The feature table containing sequence abundances per sample.'),
        'seq': ('Sequences (e.g. ASVs or representative OTUs) corresponding to '
                'the abundance table given.')
    },

    parameter_descriptions={
        'threads': 'Number of threads/processes to use during workflow.',
        'hsp_method': 'Which hidden-state prediction method to use.',
        'placement_tool': ('Placement tool to use when placing sequences into '
                           'reference tree. EPA-ng is the default, but SEPP '
                           'is less memory intensive.'),
        'min_align': ('Proportion of the total length of an input query '
                      'sequence that must align with reference sequences. '
                      'Any sequences with lengths below this value after '
                      'making an alignment with reference sequences will '
                      'be excluded from the placement and all subsequent '
                      'steps.'),
        'max_nsti': ('Max nearest-sequenced taxon index for an input ASV to '
                     'be output.'),
        'skip_minpath': ('Do not run MinPath to identify which pathways are '
                         'present as a first pass (on by default).'),
        'edge_exponent': ('Setting for maximum parisomony hidden-state '
                          'prediction. Specifies weighting transition costs '
                          'by the inverse length of edge lengths. If 0, then '
                          'edge lengths do not influence predictions. Must be '
                          'a non-negative real-valued number.'),
        'no_gap_fill': ('Do not perform gap filling before predicting '
                        'pathway abundances (gap filling is on otherwise by '
                        'default).'),
        'skip_norm': ('Skip normalizing sequence abundances by predicted '
                      'marker gene copy numbers (typically 16S rRNA '
                      'genes). The normalization step will be performed '
                      'automatically unless this option is specified.'),
        'highly_verbose': ('Print all commands being written as well as all '
                           'standard output of wrapped tools. This can be '
                           'especially useful for debugging. Note that this '
                           'option requires that the --verbose option is also '
                           'set (which is an internal QIIME 2 option that '
                           'indicates that STDOUT and STDERR should be printed '
                           'out).')
        },

    output_descriptions={'ko_metagenome': 'Predicted metagenome for KEGG orthologs',
                         'ec_metagenome': 'Predicted metagenome for EC numbers',
                         'pathway_abundance': 'Predicted MetaCyc pathway abundances'},

    name='Default 16S PICRUSt2 Pipeline',

    description=("QIIME 2 plugin for default 16S PICRUSt2 pipeline"),

    citations=[citations['Douglas2020NatureBiotech']]
)


plugin.methods.register_function(
    function=q2_picrust2.custom_tree_pipeline,

    inputs={'table': FeatureTable[Frequency],
            'tree': Phylogeny[Rooted]},

    parameters={'threads': Int % Range(1, None),
                'hsp_method': Str % Choices(HSP_METHODS),
                'max_nsti': Float % Range(0.0, None),
                'edge_exponent': Float % Range(0.0, None),
                'skip_minpath': Bool,
                'no_gap_fill': Bool,
                'skip_norm': Bool,
                'highly_verbose': Bool},

    outputs=[
       ('ko_metagenome', FeatureTable[Frequency]),
       ('ec_metagenome', FeatureTable[Frequency]),
       ('pathway_abundance', FeatureTable[Frequency])
    ],

    input_descriptions={
        'table': ('The feature table containing sequence abundances per sample.'),
        'tree': ('Tree of study ASVs placed into reference phylogeny.')
    },

    parameter_descriptions={
        'threads': 'Number of threads/processes to use during workflow.',
        'hsp_method': 'Which hidden-state prediction method to use.',
        'max_nsti': ('Max nearest-sequenced taxon index for an input ASV to '
                     'be output.'),
        'skip_minpath': ('Do not run MinPath to identify which pathways are '
                         'present as a first pass (on by default).'),
        'edge_exponent': ('Setting for maximum parisomony hidden-state '
                          'prediction. Specifies weighting transition costs '
                          'by the inverse length of edge lengths. If 0, then '
                          'edge lengths do not influence predictions. Must be '
                          'a non-negative real-valued number.'),
        'no_gap_fill': ('Do not perform gap filling before predicting '
                        'pathway abundances (gap filling is on otherwise by '
                        'default).'),
        'skip_norm': ('Skip normalizing sequence abundances by predicted '
                      'marker gene copy numbers (typically 16S rRNA '
                      'genes). The normalization step will be performed '
                      'automatically unless this option is specified.'),
        'highly_verbose': ('Print all commands being written as well as all '
                           'standard output of wrapped tools. This can be '
                           'especially useful for debugging. Note that this '
                           'option requires that the --verbose option is also '
                           'set (which is an internal QIIME 2 option that '
                           'indicates that STDOUT and STDERR should be printed '
                           'out).')
        },

    output_descriptions={'ko_metagenome': 'Predicted metagenome for KEGG orthologs',
                         'ec_metagenome': 'Predicted metagenome for E.C. numbers',
                         'pathway_abundance': 'Predicted MetaCyc pathway abundances'},

    name='16S PICRUSt2 pipeline with custom tree',

    description=("QIIME 2 plugin for running PICRUSt2 pipeline based on a " +
                 "tree from a different pipeline. This was written to be " +
                 "used with the output of SEPP (q2-fragment-insertion) as a " +
                 "starting point."),

    citations=[citations['Douglas2020NatureBiotech']]
)
