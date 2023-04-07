import skbio
import biom
from os import path
import sys
from tempfile import TemporaryDirectory
from picrust2.util import system_call_check


def custom_tree_pipeline(table: biom.Table,
                         tree: skbio.TreeNode,
                         threads: int = 1,
                         hsp_method: str = "mp",
                         max_nsti: float = 2.0,
                         edge_exponent: float = 0.5,
                         skip_minpath: bool = False,
                         no_gap_fill: bool = False,
                         skip_norm: bool = False,
                         highly_verbose: bool = False) -> (biom.Table,
                                                           biom.Table,
                                                           biom.Table):

    # Run pipeline in temporary directory so that files are not saved locally.
    with TemporaryDirectory() as temp_dir:

        # Need to write out BIOM table and newick tree to be used in pipeline.

        # Write out biom table:
        biom_infile = path.join(temp_dir, "intable.biom")
        with biom.util.biom_open(biom_infile, 'w') as out_biom:
            table.to_hdf5(h5grp=out_biom,
                          generated_by="PICRUSt2 QIIME 2 Plugin")

        # Write out newick tree.
        newick_infile = path.join(temp_dir, "placed_seqs.tre")
        tree.write(newick_infile, format="newick")

        picrust2_out = path.join(temp_dir, "picrust2_out")

        print("Running the below commands:", file=sys.stderr)

        # Run hidden-state prediction step (on 16S, EC, and KO tables
        # separately.
        hsp_out_16S = path.join(picrust2_out, "16S_predicted.tsv.gz")
        hsp_out_16S_cmd = "hsp.py -i 16S " + \
                          " -t " + newick_infile + \
                          " -p 1 " + \
                          " -n " + \
                          " -o " + hsp_out_16S + \
                          " -m " + hsp_method + \
                          " -e " + str(edge_exponent)

        hsp_out_EC = path.join(picrust2_out, "EC_predicted.tsv.gz")
        hsp_out_EC_cmd = "hsp.py -i EC " + \
                          " -t " + newick_infile + \
                          " -p " + str(threads) + \
                          " -n " + \
                          " -o " + hsp_out_EC + \
                          " -m " + hsp_method + \
                          " -e " + str(edge_exponent)

        hsp_out_KO = path.join(picrust2_out, "KO_predicted.tsv.gz")
        hsp_out_KO_cmd = "hsp.py -i KO " + \
                          " -t " + newick_infile + \
                          " -p " + str(threads) + \
                          " -n " + \
                          " -o " + hsp_out_KO + \
                          " -m " + hsp_method + \
                          " -e " + str(edge_exponent)

        if highly_verbose:
            hsp_out_16S_cmd += " --verbose"
            hsp_out_EC_cmd += " --verbose"
            hsp_out_KO_cmd += " --verbose"

        if not skip_norm:
            system_call_check(hsp_out_16S_cmd,
                              print_command=True,
                              print_stdout=highly_verbose,
                              print_stderr=True)

        system_call_check(hsp_out_EC_cmd,
                          print_command=True,
                          print_stdout=highly_verbose,
                          print_stderr=True)

        system_call_check(hsp_out_KO_cmd,
                          print_command=True,
                          print_stdout=highly_verbose,
                          print_stderr=True)

        # Run metagenome pipeline step.
        EC_metagenome_out = path.join(picrust2_out, "EC_metagenome_out")
        KO_metagenome_out = path.join(picrust2_out, "KO_metagenome_out")

        EC_metagenome_cmd = "metagenome_pipeline.py -i " + biom_infile + \
                            " -f " + hsp_out_EC + \
                            " -o " + EC_metagenome_out + \
                            " --max_nsti " + str(max_nsti)

        KO_metagenome_cmd = "metagenome_pipeline.py -i " + biom_infile + \
                            " -f " + hsp_out_KO + \
                            " -o " + KO_metagenome_out + \
                            " --max_nsti " + str(max_nsti)

        if skip_norm:
            EC_metagenome_cmd += " --skip_norm"
            KO_metagenome_cmd += " --skip_norm"
        else:
            EC_metagenome_cmd += " -m " + hsp_out_16S
            KO_metagenome_cmd += " -m " + hsp_out_16S

        system_call_check(EC_metagenome_cmd, print_command=True,
                          print_stdout=highly_verbose,
                          print_stderr=True)
        system_call_check(KO_metagenome_cmd, print_command=True,
                          print_stdout=highly_verbose,
                          print_stderr=True)

        EC_out = path.join(EC_metagenome_out, "pred_metagenome_unstrat.tsv.gz")
        KO_out = path.join(KO_metagenome_out, "pred_metagenome_unstrat.tsv.gz")

        # Run pathway inference step.
        pathways_out = path.join(picrust2_out, "pathways_out")
        pathabun_out = path.join(pathways_out, "path_abun_unstrat.tsv.gz")

        pathway_pipeline_cmd = "pathway_pipeline.py -i " + EC_out + \
                               " -o " + pathways_out + \
                               " -p " + str(threads)

        if skip_minpath:
            pathway_pipeline_cmd += " --skip_minpath"

        if no_gap_fill:
            pathway_pipeline_cmd += " --no_gap_fill"

        if highly_verbose:
            pathway_pipeline_cmd += " --verbose"

        system_call_check(pathway_pipeline_cmd, print_command=True,
                          print_stdout=highly_verbose,
                          print_stderr=True)

        # Read in output unstratified metagenome tables and return as BIOM
        # objects.
        ko_biom = biom.load_table(KO_out)
        ec_biom = biom.load_table(EC_out)
        pathabun_biom = biom.load_table(pathabun_out)

        return ko_biom, ec_biom, pathabun_biom
