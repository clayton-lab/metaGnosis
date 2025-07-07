#!/usr/bin/env bash
set -o pipefail

python_version="$(python -V | awk '{print $2}' | cut -d '.' -f 1-2)"
dram_basedir="${CONDA_PREFIX}/lib/python${python_version}/site-packages/mag_annotator"

# DRAM code is updated based on this issue: https://github.com/WrightonLabCSU/DRAM/issues/408
#TODO: Add other changes from that thread and test DRAM features

sed -i \
    -e '1043,1046d' \
    -e '1029i\
    else:\
        logger.warning(\
             "No KEGG source provided so distillation will be of limited use."\
        )' \
    -e '1608,1611c\
                              custom_hmm_cutoffs_loc=(), use_uniref=False, use_camper=False, use_vogdb=False,\
                              kofam_use_dbcan2_thresholds=False, rename_genes=True,\
                              keep_tmp_dir=True, low_mem_mode=False, threads=10, verbose=True, log_file_path: str = None, config_loc: str = None,):' \
    -e '1615,1617c\
                          custom_hmm_cutoffs_loc, use_uniref, use_camper,\
                          use_vogdb, kofam_use_dbcan2_thresholds,\
                          rename_genes, keep_tmp_dir, low_mem_mode, threads, verbose, log_file_path, config_loc,)' \
    "${dram_basedir}/annotate_bins.py"

sed -i '125,127d' "${CONDA_PREFIX}/bin/DRAM-setup.py"

# This probably isn't needed
#sed -i '250a\
#        gene_ko_link_loc=None,' \
#    "${dram_basedir}/database_handler.py"

sed -i '189,192c\
    virsorter_genes_copy = virsorter_genes.copy()\
    virsorter_genes_copy["start_position"] = virsorter_genes_copy["start_position"].astype(int)\
    virsorter_genes_copy["end_position"] = virsorter_genes_copy["end_position"].astype(int)\
    virsorter_genes = virsorter_genes_copy.sort_values("start_position")\
    virsorter_gene_number = 0' \
    "${dram_basedir}/annotate_vgfs.py"

sed -i '43,44s/kegg/uniref/' "${dram_basedir}/database_setup.py"

sed -i '503c\
        "kegg": {"gene_ko_link_loc": gene_ko_link_loc, "download_date": kegg_download_date},' \
    "${dram_basedir}/database_processing.py"
