import biom
from os import path
import pandas as pd
import skbio
from tempfile import TemporaryDirectory
import picrust2.pipeline
from picrust2.default import (default_ref_dir_bac, default_ref_dir_arc,
                                default_tables_bac, default_tables_arc,
                                default_regroup_map, default_pathway_map)

def full_pipeline(table: biom.Table,
                  seq: pd.Series,
                  threads: int = 1,
                  hsp_method: str = "mp",
                  placement_tool: str = "epa-ng",
                  min_align: float = 0.8,
                  chunk_size: int = 5000,
                  max_nsti: float = 2.0,
                  edge_exponent: float = 0.5,
                  skip_minpath: bool = False,
                  no_gap_fill: bool = False,
                  skip_norm: bool = False,
                  highly_verbose: bool = False) -> (biom.Table,
                                                    biom.Table,
                                                    biom.Table,
                                                    skbio.TreeNode,
                                                    skbio.TreeNode):

    # Write out BIOM table and FASTA to be used in pipeline.
    with TemporaryDirectory() as temp_dir:

        # Write out BIOM table:
        biom_infile = path.join(temp_dir, "intable.biom")
        with biom.util.biom_open(biom_infile, 'w') as out_biom:
            table.to_hdf5(h5grp=out_biom,
                          generated_by="PICRUSt2 QIIME 2 Plugin")

        # Write out Pandas series as FASTA:
        seq_outfile = path.join(temp_dir, "seqs.fna")

        with open(seq_outfile, "w") as outfile_fh:
            for seqname, sequence in seq.items():
                print(">" + str(seqname) + "\n" + str(sequence),
                      file=outfile_fh)

        picrust2_out = path.join(temp_dir, "picrust2_out")

        func_outputs, pathway_outputs = picrust2.pipeline.full_pipeline_split(study_fasta=seq_outfile,
                                                                        input_table=biom_infile,
                                                                        output_folder=picrust2_out,
                                                                        processes=threads,
                                                                        placement_tool=placement_tool,
                                                                        ref_dir1=default_ref_dir_bac,
                                                                        ref_dir2=default_ref_dir_arc,
                                                                        in_traits="EC,KO",
                                                                        custom_trait_tables_ref1=None,
                                                                        custom_trait_tables_ref2=None,
                                                                        marker_gene_table_ref1=default_tables_bac["16S"],
                                                                        marker_gene_table_ref2=default_tables_arc["16S"],
                                                                        pathway_map=default_pathway_map,
                                                                        rxn_func="EC",
                                                                        no_pathways=False,
                                                                        regroup_map=default_regroup_map,
                                                                        no_regroup=False,
                                                                        stratified=False,
                                                                        max_nsti=max_nsti,
                                                                        chunk_size=chunk_size,
                                                                        min_reads=1,
                                                                        min_samples=1,
                                                                        hsp_method=hsp_method,
                                                                        edge_exponent=edge_exponent,
                                                                        min_align=min_align,
                                                                        skip_minpath=skip_minpath,
                                                                        no_gap_fill=no_gap_fill,
                                                                        coverage=False,
                                                                        per_sequence_contrib=False,
                                                                        wide_table=False,
                                                                        skip_norm=skip_norm,
                                                                        remove_intermediate=False,
                                                                        verbose=highly_verbose)

        # Convert the returned unstratified tables to BIOM tables.
        # Note that the 0-index in the func table returned objects corresponds
        # to the path to the unstratified table.
        ko_biom = biom.load_table(func_outputs["KO"][0])
        ec_biom = biom.load_table(func_outputs["EC"][0])
        pathabun_biom = biom.load_table(pathway_outputs["unstrat_abun"])
        bac_tree = skbio.TreeNode.read(path.join(picrust2_out, "bac_reduced.tre"))
        arc_tree = skbio.TreeNode.read(path.join(picrust2_out, "arc_reduced.tre"))
        combined_nsti = pd.read_csv(path.join(picrust2_out, "combined_marker_predicted_and_nsti.tsv.gz"), sep="\t", index_col="sequence")

        return ko_biom, ec_biom, pathabun_biom, bac_tree, arc_tree, combined_nsti
