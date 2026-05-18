rule colabfold:
    input:
        primary_fasta=str(RESULTS / "stage2" / "proteinmpnn" / "3hpj_clean.fasta"),
        control_fasta=str(RESULTS / "stage2" / "proteinmpnn" / "2bnr_clean.fasta"),
    output:
        primary_pdb=str(RESULTS / "stage2" / "colabfold" / "3hpj_clean" / "ranked_0.pdb"),
        primary_scores=str(RESULTS / "stage2" / "colabfold" / "3hpj_clean" / "scores.json"),
        control_pdb=str(RESULTS / "stage2" / "colabfold" / "2bnr_clean" / "ranked_0.pdb"),
        control_scores=str(RESULTS / "stage2" / "colabfold" / "2bnr_clean" / "scores.json"),
    params:
        mock=MOCK,
        fixture_pdb=str(FIXTURES / "colabfold" / "sample.pdb"),
        fixture_scores=str(FIXTURES / "colabfold" / "sample_scores.json"),
        fasta_dir=str(RESULTS / "stage2" / "proteinmpnn"),
        out_dir=str(RESULTS / "stage2" / "colabfold"),
    shell:
        r"""
        mkdir -p $(dirname {output.primary_pdb}) $(dirname {output.control_pdb})
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture_pdb} {output.primary_pdb}
            cp {params.fixture_scores} {output.primary_scores}
            cp {params.fixture_pdb} {output.control_pdb}
            cp {params.fixture_scores} {output.control_scores}
        else
            python workflow/scripts/run_colabfold.py \
                --fasta-dir {params.fasta_dir} \
                --out-dir {params.out_dir}
        fi
        """
