# Cycle 03 -- Stage 2 sub-run A (truncated, 3-chain designs).
#
# Spec: specs/stage2_proteinmpnn_af2.md + PROJECT_STATE 12.j, 14, 16.
#
# Parallel to rule `stage2_designs` (02b), which handles the cycle-02 full
# 4-chain hero workflow and stays untouched. Sub-run A designs are 3-chain
# (A=binder, B=HLA[1:180], C=peptide); this rule runs run_stage2.py with
# --target-layout truncated, which forks the splice (splice_binder_subrun_a ->
# A=HLA, B=peptide, C=binder), designs ProteinMPNN chain C, builds the AF2 FASTA
# in HLA:peptide:binder order, and decomposes metrics via
# LAYOUT_CHAINS["truncated"]. beta2m is never reintroduced in silico.
#
# This is an explicit-target rule (not part of `rule all`): the default cycle-03
# build still drives the 4-chain stage2_designs, while the truncated funnel is
# requested directly, e.g.
#   snakemake results/cycle_03/stage2/subrun_a/stage2_summary.json --config cycle=03 mock=false -j1
#
# Mock mode uses the dedicated 3-chain truncated fixtures
# (tests/fixtures/stage2/designs_truncated): the stage1_subrun_a mock synthesizes
# 4-chain designs, which the truncated splice deliberately rejects.

rule stage2_subrun_a:
    input:
        subrun_summary=str(RESULTS / "stage1" / "subrun_a" / "subrun_summary.json"),
        subrun_designs=str(RESULTS / "stage1" / "subrun_a" / "designs.jsonl"),
        target_yaml=str(RESULTS / "stage0" / "target.yaml"),
        truncated="data/targets/3hpj_baker_truncated.pdb",
        mpnn_cfg=MPNN_DESIGNS_CONFIG,
        af2_cfg="configs/af2_stage2.yaml",
        seeds="configs/seeds.yaml",
    output:
        sequences_jsonl=str(
            RESULTS / "stage2" / "subrun_a" / "proteinmpnn_designs" / "sequences.jsonl"
        ),
        spliced_marker=str(
            RESULTS / "stage2" / "subrun_a" / "proteinmpnn_designs" / "spliced" / ".done"
        ),
        metrics_jsonl=str(RESULTS / "stage2" / "subrun_a" / "af2_designs" / "metrics.jsonl"),
        summary=str(RESULTS / "stage2" / "subrun_a" / "stage2_summary.json"),
    params:
        mock=MOCK,
        cycle=CYCLE,
        out_dir=str(RESULTS / "stage2" / "subrun_a"),
        mock_summary="tests/fixtures/stage2/designs_truncated/subrun_summary.json",
        mock_manifest="tests/fixtures/stage2/designs_truncated/mock_truncated_target.yaml",
    shell:
        r"""
        mkdir -p {params.out_dir}
        if [ "{params.mock}" = "True" ]; then
            # TODO(Trap #42): this mock branch reads the dedicated 3-chain
            # truncated fixtures (params.mock_summary/mock_manifest) instead of the
            # declared input.subrun_summary, because the stage1_subrun_a mock still
            # synthesizes 4-chain designs (which the truncated splice rejects).
            # Re-couple to input.subrun_summary once stage1_subrun_a mock produces
            # truncated 3-chain designs. The real branch below already reads the
            # declared inputs correctly.
            python scripts/run_stage2.py \
                --mock --cycle {params.cycle} \
                --target-layout truncated \
                --stage1-summary {params.mock_summary} \
                --target-manifest {params.mock_manifest} \
                --out-dir {params.out_dir}
        else
            python scripts/run_stage2.py \
                --cycle {params.cycle} \
                --target-layout truncated \
                --stage1-summary {input.subrun_summary} \
                --target-manifest {input.target_yaml} \
                --cleaned-pmhc {input.truncated} \
                --mpnn-config {input.mpnn_cfg} \
                --af2-config {input.af2_cfg} \
                --seeds-config {input.seeds} \
                --out-dir {params.out_dir}
        fi
        """
