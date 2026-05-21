# Cycle 03 -- Stage 1 sub-run B: design_2079 hero partial diffusion.
#
# Spec: specs/stage1_rfdiffusion.md ("Cycle 03: partial-diffusion subruns").
# Partial-diffuses num_designs binder backbones from the cycle-02 hero seed
# complex (data/seeds/design_2079_binder.pdb, built pod-side by
# setup_cycle03_inputs.py). Active only when stage1_mode=partial (cycle >= 3).

rule stage1_subrun_b:
    input:
        target_yaml=str(RESULTS / "stage0" / "target.yaml"),
        cfg="configs/rfdiffusion_subrun_b.yaml",
        seeds="configs/seeds.yaml",
    output:
        summary=str(RESULTS / "stage1" / "subrun_b" / "subrun_summary.json"),
        designs_jsonl=str(RESULTS / "stage1" / "subrun_b" / "designs.jsonl"),
    params:
        mock=MOCK,
        cycle=CYCLE,
        out_dir=str(RESULTS / "stage1" / "subrun_b"),
    shell:
        r"""
        mkdir -p {params.out_dir}
        if [ "{params.mock}" = "True" ]; then
            python scripts/run_stage1_subrun.py \
                --mock --subrun b --cycle {params.cycle} \
                --out-dir {params.out_dir}
        else
            python scripts/run_stage1_subrun.py \
                --subrun b --cycle {params.cycle} \
                --config {input.cfg} \
                --seeds-yaml {input.seeds} \
                --target-manifest {input.target_yaml} \
                --out-dir {params.out_dir}
        fi
        """
