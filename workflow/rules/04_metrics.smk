rule metrics:
    input:
        scores=str(RESULTS / "stage2" / "colabfold" / "design_00000" / "scores.json"),
        thresholds="configs/thresholds.yaml",
    output:
        metrics=str(RESULTS / "stage2" / "metrics.json"),
        survivors=str(RESULTS / "stage2" / "survivors.json"),
    params:
        mock=MOCK,
        fixture=str(FIXTURES / "metrics" / "sample_metrics.json"),
    shell:
        r"""
        mkdir -p $(dirname {output.metrics})
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture} {output.metrics}
            printf '{{"mock": true, "survivors": ["design_00000"]}}\n' > {output.survivors}
        else
            python workflow/scripts/compute_metrics.py --config {input.scores} --out {output.metrics}
        fi
        """
