import qiime2
import biom
from os import path
import pandas as pd
from tempfile import TemporaryDirectory
from q2_types.feature_table import FeatureTable, Frequency
import subprocess
import sys
import picrust2.pipeline
from picrust2.default import (default_ref_dir, default_tables, default_map,
                              default_regroup_map, default_pathway_map)

def full_pipeline(table: biom.Table,
                  seq: pd.Series,
                  threads: int = 1,
                  hsp_method: str = "mp",
                  max_nsti: float = 2.0) -> (biom.Table,
                                             biom.Table,
                                             biom.Table):

    # Write out BIOM table and FASTA to be used in pipeline.
    with TemporaryDirectory() as temp_dir:

        # Write out BIOM table:
        biom_infile = path.join(temp_dir, "intable.biom")
        with biom.util.biom_open(biom_infile, 'w') as out_biom:  
            table.to_hdf5(h5grp=out_biom,
                          generated_by="PICRUSt2 QIIME2 Plugin")

        # Write out Pandas series as FASTA:
        seq_outfile = path.join(temp_dir, "seqs.fna")

        with open(seq_outfile, "w") as outfile_fh:
            for seqname, sequence in seq.iteritems():
                print(">" + str(seqname) + "\n" + str(sequence),
                      file=outfile_fh)

        picrust2_out = path.join(temp_dir, "picrust2_out")

        func_outputs, pathway_outputs = picrust2.pipeline.full_pipeline(study_fasta=seq_outfile,
                                                                        input_table=biom_infile,
                                                                        output_folder=picrust2_out,
                                                                        processes=threads,
                                                                        ref_dir=default_ref_dir,
                                                                        in_traits="EC,KO",
                                                                        custom_trait_tables=None,
                                                                        marker_gene_table=default_tables["16S"],
                                                                        pathway_map=default_pathway_map,
                                                                        rxn_func="EC",
                                                                        no_pathways=False,
                                                                        regroup_map=default_regroup_map,
                                                                        no_regroup=False,
                                                                        metagenome_contrib=True,
                                                                        stratified=False,
                                                                        max_nsti=max_nsti,
                                                                        min_reads=1,
                                                                        min_samples=1,
                                                                        hsp_method=hsp_method,
                                                                        skip_nsti=False,
                                                                        no_gap_fill=False,
                                                                        skip_minpath=False,
                                                                        coverage=False,
                                                                        per_sequence_contrib=False,
                                                                        remove_intermediate=False,
                                                                        verbose=True)

        # Convert the returned unstratified tables to BIOM tables.
        # Note that the 0-index in the func table returned objects corresponds
        # to the path to the unstratified table.
        ko_biom = biom.load_table(func_outputs["KO"][0])
        ec_biom = biom.load_table(func_outputs["EC"][0])
        pathabun_biom = biom.load_table(pathway_outputs["unstrat_abun"])

        return ko_biom, ec_biom, pathabun_biom
