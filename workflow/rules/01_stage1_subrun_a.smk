# Cycle 03 -- Stage 1 sub-run A: BAKER scaffold library partial diffusion.
#
# Spec: specs/stage1_rfdiffusion.md ("Cycle 03: partial-diffusion subruns").
# Aligns each BAKER scaffold onto our reference (align_scaffolds.py) then
# partial-diffuses one binder backbone per scaffold (partial_diffuse.py).
# Active only when stage1_mode=partial (cycle >= 3); see Snakefile ruleorder.

rule stage1_subrun_a:
    input:
        target_yaml=str(RESULTS / "stage0" / "target.yaml"),
        cfg="configs/rfdiffusion_subrun_a.yaml",
        seeds="configs/seeds.yaml",
    output:
        summary=str(RESULTS / "stage1" / "subrun_a" / "subrun_summary.json"),
        designs_jsonl=str(RESULTS / "stage1" / "subrun_a" / "designs.jsonl"),
    params:
        mock=MOCK,
        cycle=CYCLE,
        out_dir=str(RESULTS / "stage1" / "subrun_a"),
    shell:
        r"""
        mkdir -p {params.out_dir}
        if [ "{params.mock}" = "True" ]; then
            python scripts/run_stage1_subrun.py \
                --mock --subrun a --cycle {params.cycle} \
                --out-dir {params.out_dir}
        else
            python scripts/run_stage1_subrun.py \
                --subrun a --cycle {params.cycle} \
                --config {input.cfg} \
                --seeds-yaml {input.seeds} \
                --target-manifest {input.target_yaml} \
                --out-dir {params.out_dir}
        fi
        """
