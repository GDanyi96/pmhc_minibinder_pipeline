rule colabfold:
    input:
        fasta=str(RESULTS / "stage2" / "proteinmpnn" / "design_00000.fasta"),
    output:
        pdb=str(RESULTS / "stage2" / "colabfold" / "design_00000" / "ranked_0.pdb"),
        scores=str(RESULTS / "stage2" / "colabfold" / "design_00000" / "scores.json"),
    params:
        mock=MOCK,
        fixture_pdb=str(FIXTURES / "colabfold" / "sample.pdb"),
        fixture_scores=str(FIXTURES / "colabfold" / "sample_scores.json"),
    shell:
        r"""
        mkdir -p $(dirname {output.pdb})
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture_pdb} {output.pdb}
            cp {params.fixture_scores} {output.scores}
        else
            python workflow/scripts/run_colabfold.py --config {input.fasta} --out {output.scores}
        fi
        """
