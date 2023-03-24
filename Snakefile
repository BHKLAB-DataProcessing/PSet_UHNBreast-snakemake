from os import path
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"],
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]
filename = config["filename"]
rna_tool = config["rna_tool"]
rna_ref = config["rna_ref"]
is_filtered = config["filtered"]
filtered = filtered = 'filtered' if is_filtered is not None and is_filtered == 'True' else ''

rna_tool_dir = rna_tool.replace('-', '_')
rnaseq_dir = path.join(prefix, "processed",
                       rna_tool_dir, rna_tool_dir + '_' + rna_ref)
rna_ref_file = rna_ref.replace('_', '.') + '.annotation.RData'

rule get_pset:
    input:
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv"),
        S3.remote(prefix + "download/UHN_recomputed.RData"),
        S3.remote(prefix + "download/" + rna_tool_dir + '.tar.gz')
    output:
        S3.remote(prefix + filename)
    shell:
        """
        Rscript {prefix}scripts/UHNBreast_2019.R {prefix} {rna_tool} {rna_ref} {filtered}
        """

rule download_annotation:
    output:
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv"),
        S3.remote(prefix + 'download/bc_cellines_neel_subtypes.csv'),
        S3.remote(prefix + 'download/' + rna_ref_file),
        S3.remote(prefix + 'download/uhn_metadata_new.csv')
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
        S3.remote(prefix + "download/UHN_recomputed.RData"),
        S3.remote(prefix + "download/" + rna_tool_dir + '.tar.gz')
    shell:
        """
        Rscript {prefix}scripts/download_data.R {prefix} {rna_tool_dir}.tar.gz
        """
