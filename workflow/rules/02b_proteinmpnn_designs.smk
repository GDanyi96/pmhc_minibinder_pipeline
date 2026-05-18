# Stage 2 designs (real Stage 1 outputs) -- ProteinMPNN sequence design +
# splice + AF2 fan-in + halt gate. Single rule runs the orchestrator end
# to end; downstream 03b is a marker-only DAG anchor.
#
# Spec: specs/stage2_proteinmpnn_af2.md
#
# Distinct from rule `proteinmpnn` (in 02_proteinmpnn.smk) which handles
# the P1-N2 controls panel. The two pipelines share the validated MPNN +
# ColabFold + compute_metrics drivers but have independent halt gates.

rule stage2_designs:
    input:
        stage1_summary=str(RESULTS / "stage1" / "rfdiffusion" / "stage1_summary.json"),
        target_yaml=str(RESULTS / "stage0" / "target.yaml"),
        mpnn_cfg="configs/proteinmpnn_default.yaml",
        af2_cfg="configs/af2_stage2.yaml",
        seeds="configs/seeds.yaml",
    output:
        sequences_jsonl=str(RESULTS / "stage2" / "proteinmpnn_designs" / "sequences.jsonl"),
        spliced_marker=str(RESULTS / "stage2" / "proteinmpnn_designs" / "spliced" / ".done"),
        metrics_jsonl=str(RESULTS / "stage2" / "af2_designs" / "metrics.jsonl"),
        summary=str(RESULTS / "stage2" / "stage2_summary.json"),
    params:
        mock=MOCK,
        cycle=CYCLE,
        out_dir=str(RESULTS / "stage2"),
    shell:
        r"""
        mkdir -p {params.out_dir}
        if [ "{params.mock}" = "True" ]; then
            python scripts/run_stage2.py \
                --mock --cycle {params.cycle} \
                --out-dir {params.out_dir}
        else
            python scripts/run_stage2.py \
                --cycle {params.cycle} \
                --stage1-summary {input.stage1_summary} \
                --target-manifest {input.target_yaml} \
                --mpnn-config {input.mpnn_cfg} \
                --af2-config {input.af2_cfg} \
                --seeds-config {input.seeds} \
                --out-dir {params.out_dir}
        fi
        """
