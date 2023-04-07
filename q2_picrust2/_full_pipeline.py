import biom
from os import path
import pandas as pd
from tempfile import TemporaryDirectory
import picrust2.pipeline
from picrust2.default import (default_ref_dir, default_tables,
                              default_regroup_map, default_pathway_map)


def full_pipeline(table: biom.Table,
                  seq: pd.Series,
                  threads: int = 1,
                  hsp_method: str = "mp",
                  placement_tool: str = "epa-ng",
                  min_align: float = 0.8,
                  max_nsti: float = 2.0,
                  edge_exponent: float = 0.5,
                  skip_minpath: bool = False,
                  no_gap_fill: bool = False,
                  skip_norm: bool = False,
                  highly_verbose: bool = False) -> (biom.Table,
                                                    biom.Table,
                                                    biom.Table):

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

        func_outputs, pathway_outputs = picrust2.pipeline.full_pipeline(study_fasta=seq_outfile,
                                                                        input_table=biom_infile,
                                                                        output_folder=picrust2_out,
                                                                        processes=threads,
                                                                        placement_tool=placement_tool,
                                                                        ref_dir=default_ref_dir,
                                                                        in_traits="EC,KO",
                                                                        custom_trait_tables=None,
                                                                        marker_gene_table=default_tables["16S"],
                                                                        pathway_map=default_pathway_map,
                                                                        rxn_func="EC",
                                                                        no_pathways=False,
                                                                        regroup_map=default_regroup_map,
                                                                        no_regroup=False,
                                                                        stratified=False,
                                                                        max_nsti=max_nsti,
                                                                        min_reads=1,
                                                                        min_samples=1,
                                                                        hsp_method=hsp_method,
                                                                        edge_exponent=edge_exponent,
                                                                        min_align=min_align,
                                                                        skip_nsti=False,
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

        return ko_biom, ec_biom, pathabun_biom
