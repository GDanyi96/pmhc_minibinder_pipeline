rule prep_target:
    input:
        cfg="configs/target_wt1_a0201.yaml",
    output:
        target_yaml=str(RESULTS / "stage0" / "target.yaml"),
        cleaned_pdb=str(RESULTS / "stage0" / "target_clean.pdb"),
    params:
        mock=MOCK,
        fixture=str(FIXTURES / "target_3hpj_clean.pdb"),
    shell:
        r"""
        mkdir -p $(dirname {output.target_yaml})
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture} {output.cleaned_pdb}
            printf 'mock: true\ntarget_id: wt1_a0201\n' > {output.target_yaml}
        else
            python workflow/scripts/prep_target.py --config {input.cfg} --out {output.target_yaml}
        fi
        """
