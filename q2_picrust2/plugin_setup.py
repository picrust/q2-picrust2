from qiime2.plugin import (Plugin, Str, Properties, Choices, Int, Bool, Range,
                           Float, Set, Visualization, Metadata, MetadataColumn,
                           Categorical, Numeric, Citations)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Sequence
from q2_types.sample_data import SampleData
from q2_types.tree import Phylogeny, Rooted
import q2_picrust2

citations = Citations.load('citations.bib', package='q2_picrust2')

HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']

plugin = Plugin(
    name='picrust2',
    version="0.0.3",
    website='https://github.com/gavinmdouglas/q2-picrust2',
    package='q2_picrust2',
    description=('This QIIME 2 plugin wraps the default 16S PICRUSt2 pipeline to run '
                 'metagenome inference based on marker gene data. Currently '
                 'only unstratified output is supported.'),
    short_description='Predicts gene families and pathways from 16S sequences.',
)

plugin.methods.register_function(
    function=q2_picrust2.full_pipeline,

    inputs={'table': FeatureTable[Frequency],
            'seq': FeatureData[Sequence]},
             
    parameters={'threads': Int % Range(1, None),
                'hsp_method': Str % Choices(HSP_METHODS),
                'max_nsti': Float % Range(0.0, None)},

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
        'max_nsti': ('Max nearest-sequenced taxon index for an input ASV to '
                     'be output.')},

    output_descriptions={'ko_metagenome': 'Predicted metagenome for KEGG orthologs',
                         'ec_metagenome': 'Predicted metagenome for E.C. numbers',
                         'pathway_abundance': 'Predicted MetaCyc pathway abundances'},

    name='Default 16S PICRUSt2 Pipeline',

    description=("QIIME2 Plugin for default 16S PICRUSt2 pipeline"),

    citations=[citations['Langille2013NatureBioTech']]
)


plugin.methods.register_function(
    function=q2_picrust2.custom_tree_pipeline,

    inputs={'table': FeatureTable[Frequency],
            'tree': Phylogeny[Rooted]},
             
    parameters={'threads': Int % Range(1, None),
                'hsp_method': Str % Choices(HSP_METHODS),
                'max_nsti': Float % Range(0.0, None)},

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
                     'be output.')},

    output_descriptions={'ko_metagenome': 'Predicted metagenome for KEGG orthologs',
                         'ec_metagenome': 'Predicted metagenome for E.C. numbers',
                         'pathway_abundance': 'Predicted MetaCyc pathway abundances'},

    name='16S PICRUSt2 pipeline with custom tree',

    description=("QIIME2 plugin for running PICRUSt2 pipeline based on a " +
                 "tree from a different pipeline. This was written to be " +
                 "used with the output of SEPP (q2-fragment-insertion) as a " +
                 "starting point."),

    citations=[citations['Langille2013NatureBioTech']]
)

