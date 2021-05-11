import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_download():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/download/data")
        expected_path = PurePosixPath(".tests/unit/download/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.dbsnp138.vcf Homo_sapiens_assembly38.known_indels.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_omni2.5.hg38.vcf.gz hapmap_3.3.hg38.vcf.gz", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.dbsnp138.vcf Homo_sapiens_assembly38.known_indels.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_omni2.5.hg38.vcf.gz hapmap_3.3.hg38.vcf.gz",
            "-F", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
