# Cycle 03 -- Stage 1 merge: combine sub-runs A + B into the canonical output.
#
# Spec: specs/stage1_rfdiffusion.md ("Cycle 03: partial-diffusion subruns").
# Produces the SAME canonical outputs as `rule rfdiffusion`
# (stage1/rfdiffusion/{stage1_summary.json, designs.jsonl}); the Snakefile's
# conditional `ruleorder` selects this rule when stage1_mode=partial (cycle
# >= 3) and `rule rfdiffusion` otherwise. The geometry halt gate runs once
# here on the merged design set.

rule stage1_merge:
    input:
        subrun_a=str(RESULTS / "stage1" / "subrun_a" / "subrun_summary.json"),
        subrun_b=str(RESULTS / "stage1" / "subrun_b" / "subrun_summary.json"),
        seeds="configs/seeds.yaml",
    output:
        summary=str(RESULTS / "stage1" / "rfdiffusion" / "stage1_summary.json"),
        designs_jsonl=str(RESULTS / "stage1" / "rfdiffusion" / "designs.jsonl"),
    params:
        mock=MOCK,
        cycle=CYCLE,
        out_dir=str(RESULTS / "stage1" / "rfdiffusion"),
        subrun_a_dir=str(RESULTS / "stage1" / "subrun_a"),
        subrun_b_dir=str(RESULTS / "stage1" / "subrun_b"),
    shell:
        r"""
        mkdir -p {params.out_dir}
        MOCK_FLAG=""
        if [ "{params.mock}" = "True" ]; then MOCK_FLAG="--mock"; fi
        python scripts/merge_stage1_subruns.py \
            $MOCK_FLAG --cycle {params.cycle} \
            --seeds-yaml {input.seeds} \
            --subrun-dir {params.subrun_a_dir} \
            --subrun-dir {params.subrun_b_dir} \
            --out-dir {params.out_dir}
        """
