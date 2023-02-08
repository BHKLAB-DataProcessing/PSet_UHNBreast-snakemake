from os import path
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(
#     access_key_id=config["key"],
#     secret_access_key=config["secret"],
#     host=config["host"],
#     stay_on_remote=False
# )

prefix = config["prefix"]
rna_tool = 'Kallisto-0.46.1'
rna_ref = 'Gencode_v33'

rna_tool_dir = rna_tool.replace('-', '_')
rnaseq_dir = path.join(prefix, "processed",
                       rna_tool_dir, rna_tool_dir + '_' + rna_ref)
rna_ref_file = rna_ref.replace('_', '.') + '.annotation.RData'

rule get_pset:
    input:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + "download/UHN_recomputed.RData",
        prefix + "download/" + rna_tool_dir + '.tar.gz'
    output:
        prefix + "UHNBreast.rds"
    shell:
        """
        Rscript {prefix}scripts/UHNBreast_2019.R {prefix} {rna_tool} {rna_ref}
        """

rule download_annotation:
    output:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + 'download/bc_cellines_neel_subtypes.csv',
        prefix + 'download/' + rna_ref_file,
        prefix + 'download/uhn_metadata_new.csv'
    shell:
        """
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/drugs_with_ids.csv' \
            -O {prefix}download/drugs_with_ids.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/cell_annotation_all.csv' \
            -O {prefix}download/cell_annotation_all.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/bc_cellines_neel_subtypes.csv' \
            -O {prefix}download/bc_cellines_neel_subtypes.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/{rna_ref_file}' \
            -O {prefix}download/{rna_ref_file}
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/uhn_metadata_new.csv' \
            -O {prefix}download/uhn_metadata_new.csv
        """

rule download_data:
    output:
        prefix + "download/UHN_recomputed.RData",
        prefix + "download/" + rna_tool_dir + '.tar.gz'
    shell:
        """
        Rscript {prefix}scripts/download_data.R {prefix} {rna_tool_dir}.tar.gz
        """
