import qiime2
import biom
from os import path
import pandas as pd
from tempfile import TemporaryDirectory
from q2_types.feature_table import FeatureTable, Frequency
import picrust2.pipeline
from picrust2.default import (default_fasta, default_tree, default_hmm,
                              default_tables, default_map, default_regroup_map,
                              default_pathway_map)

def full_pipeline(table: biom.Table,
                  seq : pd.Series,
                  threads: int = 1,
                  hsp_method: str = "mp")  -> (biom.Table,
                                               biom.Table,
                                               biom.Table,
                                               biom.Table):

    # Need to write out BIOM table and fasta to be used in pipeline.
    with TemporaryDirectory() as temp_dir:
            
        # Write out biom table:
        biom_infile = path.join(temp_dir, "intable.biom")
        with biom.util.biom_open(biom_infile, 'w') as out_biom:  
            table.to_hdf5(h5grp=out_biom, generated_by="PICRUSt2 QIIME2 Plugin")

        # Write out Pandas series as fasta:
        seq_outfile = path.join(temp_dir, "seqs.fna")

        with open(seq_outfile, "w") as outfile_fh:
            for seqname, sequence in seq.iteritems():
                print(">" + str(seqname) + "\n" + str(sequence),
                      file=outfile_fh)

        picrust2_out = path.join(temp_dir, "picrust2_out")

        func_outputs, pathway_outputs = picrust2.pipeline.full_pipeline(study_fasta=seq_outfile,
                                                                        input_table=biom_infile,
                                                                        output_folder=picrust2_out,
                                                                        threads=threads,
                                                                        ref_msa=default_fasta,
                                                                        tree=default_tree,
                                                                        hmm=default_hmm,
                                                                        in_traits="EC,KO",
                                                                        custom_trait_tables=None,
                                                                        marker_gene_table=default_tables["16S"],
                                                                        pathway_map=default_pathway_map,
                                                                        no_pathways=False,
                                                                        regroup_map=default_regroup_map,
                                                                        no_regroup=False,
                                                                        stratified=False,
                                                                        alignment_tool="hmmalign",
                                                                        max_nsti=2,
                                                                        min_reads=1,
                                                                        min_samples=1,
                                                                        hsp_method=hsp_method,
                                                                        calculate_NSTI=True,
                                                                        confidence=False,
                                                                        seed=198,
                                                                        no_gap_fill=False,
                                                                        per_sequence_contrib=False,
                                                                        no_descrip=True,
                                                                        verbose=False)

        Convert the returned unstratified tables to biom tables.
        ko_biom = biom.load_table(func_outputs["KO"])
        ec_biom = biom.load_table(func_outputs["EC"])
        pathabun_biom = biom.load_table(pathway_outputs["unstrat_abun"])
        pathcov_biom = biom.load_table(pathway_outputs["unstrat_cov"])

        return ko_biom, ec_biom, pathabun_biom, pathcov_biom
