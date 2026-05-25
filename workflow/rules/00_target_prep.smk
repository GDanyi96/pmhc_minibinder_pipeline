rule prep_target:
    input:
        primary_cfg="configs/target_wt1_a0201.yaml",
        control_cfg="configs/target_2bnr_a0201.yaml",
    output:
        target_yaml=str(RESULTS / "stage0" / "target.yaml"),
        primary_clean="data/targets/3hpj_clean.pdb",
        control_clean="data/targets/2bnr_clean.pdb",
    params:
        mock=MOCK,
        primary_fix=str(FIXTURES / "target_3hpj_clean.pdb"),
        control_fix=str(FIXTURES / "target_2bnr_clean.pdb"),
        primary_raw="data/targets/3HPJ.pdb",
        control_raw="data/targets/2BNR.pdb",
    shell:
        r"""
        mkdir -p $(dirname {output.target_yaml}) data/targets
        if [ "{params.mock}" = "True" ]; then
            cp {params.primary_fix} {output.primary_clean}
            cp {params.control_fix} {output.control_clean}
            python workflow/scripts/prep_target.py \
                --primary-config {input.primary_cfg} \
                --control-config {input.control_cfg} \
                --primary-pdb-out {output.primary_clean} \
                --control-pdb-out {output.control_clean} \
                --out {output.target_yaml} --mock
        else
            python workflow/scripts/prep_target.py \
                --primary-config {input.primary_cfg} \
                --control-config {input.control_cfg} \
                --primary-pdb-in  {params.primary_raw} \
                --control-pdb-in  {params.control_raw} \
                --primary-pdb-out {output.primary_clean} \
                --control-pdb-out {output.control_clean} \
                --out {output.target_yaml}
        fi
        """


# Cycle 03 sub-run A only: BAKER-format truncated target (chain B = HLA[1:180],
# chain C = peptide[1:9]; no beta2m, no alpha3). See prep_baker_target.py and
# specs/stage0_target_specification.md. Requested only by stage1_subrun_a, which
# exists only for cycle >= 3 (partial diffusion); cycle 01/02 never build it.
rule prep_baker_target:
    input:
        primary_clean="data/targets/3hpj_clean.pdb",
    output:
        truncated="data/targets/3hpj_baker_truncated.pdb",
    params:
        mock=MOCK,
    shell:
        r"""
        mkdir -p data/targets
        if [ "{params.mock}" = "True" ]; then
            python workflow/scripts/prep_baker_target.py --mock --out {output.truncated}
        else
            python workflow/scripts/prep_baker_target.py \
                --in {input.primary_clean} --out {output.truncated}
        fi
        """
