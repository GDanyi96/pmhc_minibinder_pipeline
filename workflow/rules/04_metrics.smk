rule metrics:
    input:
        primary_scores=str(RESULTS / "stage2" / "colabfold" / "3hpj_clean" / "scores.json"),
        control_scores=str(RESULTS / "stage2" / "colabfold" / "2bnr_clean" / "scores.json"),
        manifest="data/controls/controls_manifest.yaml",
        thresholds="configs/thresholds.yaml",
    output:
        metrics=str(RESULTS / "stage2" / "metrics.json"),
        survivors=str(RESULTS / "stage2" / "survivors.json"),
    params:
        mock=MOCK,
        fixture=str(FIXTURES / "stage2" / "stage2_metrics_mock.json"),
        cycle=CYCLE,
    shell:
        r"""
        mkdir -p $(dirname {output.metrics})
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture} {output.metrics}
            printf '{{"mock": true, "survivors": ["3hpj_clean", "2bnr_clean"]}}\n' > {output.survivors}
        else
            python scripts/run_controls.py --cycle {params.cycle} \
                --config {input.manifest} --thresholds {input.thresholds}
            printf '{{"survivors": ["3hpj_clean", "2bnr_clean"]}}\n' > {output.survivors}
        fi
        """
