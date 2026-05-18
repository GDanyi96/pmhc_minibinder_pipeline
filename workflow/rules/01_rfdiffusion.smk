# Stage 1 -- RFdiffusion backbone generation.
#
# Spec: specs/stage1_rfdiffusion.md
# Mock-mode copies tests/fixtures/stage1/ designs and a synthetic target
# manifest; real mode invokes scripts/run_stage1.py which shells out to
# /workspace/RFdiffusion/scripts/run_inference.py (native install on the
# pod -- no Docker, no DinD).

rule rfdiffusion:
    input:
        target_yaml=str(RESULTS / "stage0" / "target.yaml"),
        cfg="configs/rfdiffusion_default.yaml",
        seeds="configs/seeds.yaml",
    output:
        summary=str(RESULTS / "stage1" / "rfdiffusion" / "stage1_summary.json"),
        designs_jsonl=str(RESULTS / "stage1" / "rfdiffusion" / "designs.jsonl"),
        sample_pdb=str(RESULTS / "stage1" / "rfdiffusion" / "designs" / "design_00000.pdb"),
    params:
        mock=MOCK,
        cycle=CYCLE,
        out_dir=str(RESULTS / "stage1" / "rfdiffusion"),
        mock_manifest=str(FIXTURES / "stage1" / "mock_target.yaml"),
    shell:
        r"""
        mkdir -p $(dirname {output.summary})
        if [ "{params.mock}" = "True" ]; then
            python scripts/run_stage1.py \
                --mock --cycle {params.cycle} \
                --target-manifest {params.mock_manifest} \
                --out-dir {params.out_dir}
        else
            python scripts/run_stage1.py \
                --cycle {params.cycle} \
                --target-manifest {input.target_yaml} \
                --rfdiff-yaml {input.cfg} \
                --seeds-yaml {input.seeds} \
                --out-dir {params.out_dir}
        fi
        """
