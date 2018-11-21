import qiime2
import skbio
import biom
from os import path
import pandas as pd
from tempfile import TemporaryDirectory
from q2_types.feature_table import FeatureTable, Frequency
from picrust2.util import system_call_check

def custom_tree_pipeline(table: biom.Table,
                         tree: skbio.TreeNode,
                         threads: int = 1,
                         hsp_method: str = "mp",
                         max_nsti: float = 2.0)  -> (biom.Table,
                                                     biom.Table,
                                                     biom.Table,
                                                     biom.Table):

    # Run pipeline in temporary directory so that files are not saved locally.
    with TemporaryDirectory() as temp_dir:

        # Need to write out BIOM table and newick tree to be used in pipeline.

        # Write out biom table:
        biom_infile = path.join(temp_dir, "intable.biom")
        with biom.util.biom_open(biom_infile, 'w') as out_biom:  
            table.to_hdf5(h5grp=out_biom, generated_by="PICRUSt2 QIIME2 Plugin")

        # Write out newick tree.
        newick_infile = path.join(temp_dir, "placed_seqs.tre")
        tree.write(newick_infile, format="newick")

        picrust2_out = path.join(temp_dir, "picrust2_out")

        # Run hidden-state prediction step (on 16S, EC, and KO tables
        # separately.
        hsp_out_16S = path.join(picrust2_out, "16S_predicted")
        system_call_check("hsp.py -i 16S -t " + newick_infile + " -p 1 -n "
                          "-o " + hsp_out_16S + " -m " + hsp_method)

        hsp_out_EC = path.join(picrust2_out, "EC_predicted")
        system_call_check("hsp.py -i EC -t " + newick_infile + " -p " +
                          str(threads) + " -o " + hsp_out_EC + " -m " +
                          hsp_method)

        hsp_out_KO = path.join(picrust2_out, "KO_predicted")
        system_call_check("hsp.py -i KO -t " + newick_infile + " -p " +
                          str(threads) + " -o " + hsp_out_KO + " -m " +
                          hsp_method)

        # Run metagenome pipeline step.
        EC_metagenome_out = path.join(picrust2_out, "EC_metagenome_out")
        system_call_check("metagenome_pipeline.py -i " + biom_infile + " -m " +
                          hsp_out_16S + ".tsv -f " + hsp_out_EC + ".tsv -p " +
                          str(threads) + " -o " + EC_metagenome_out +
                          " --max_nsti " + str(max_nsti))

        KO_metagenome_out = path.join(picrust2_out, "KO_metagenome_out")
        system_call_check("metagenome_pipeline.py -i " + biom_infile + " -m " +
                          hsp_out_16S + ".tsv -f " + hsp_out_KO + ".tsv -p " +
                          str(threads) + " -o " + KO_metagenome_out +
                          " --max_nsti " + str(max_nsti))

        # Run pathway inference step.
        pathways_out = path.join(picrust2_out, "pathways_out")

        EC_out = path.join(EC_metagenome_out, "pred_metagenome_unstrat.tsv")

        system_call_check("run_minpath.py -i " + EC_out + " -o " +
                          pathways_out + " -p " + str(threads))

        # Read in output unstratified metagenome tables and return as BIOM
        # objects.
        KO_out = path.join(KO_metagenome_out, "pred_metagenome_unstrat.tsv")
        pathabun_out = path.join(pathways_out, "path_abun_unstrat.tsv")
        pathcov_out = path.join(pathways_out, "path_cov_unstrat.tsv")

        ko_biom = biom.load_table(KO_out)
        ec_biom = biom.load_table(EC_out)
        pathabun_biom = biom.load_table(pathabun_out)
        pathcov_biom = biom.load_table(pathcov_out)

        return ko_biom, ec_biom, pathabun_biom, pathcov_biom
