rule proteinmpnn:
    input:
        manifest=str(RESULTS / "stage1" / "manifest.json"),
    output:
        sample_fasta=str(RESULTS / "stage2" / "proteinmpnn" / "design_00000.fasta"),
    params:
        mock=MOCK,
        fixture=str(FIXTURES / "proteinmpnn" / "sample.fasta"),
    shell:
        r"""
        mkdir -p $(dirname {output.sample_fasta})
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture} {output.sample_fasta}
        else
            python workflow/scripts/run_proteinmpnn.py --config {input.manifest} --out {output.sample_fasta}
        fi
        """
