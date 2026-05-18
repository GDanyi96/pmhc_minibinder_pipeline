rule crosspan:
    input:
        survivors=str(RESULTS / "stage2" / "survivors.json"),
        thresholds="configs/thresholds.yaml",
    output:
        panel=str(RESULTS / "stage3" / "crosspan_panel.json"),
        spec_survivors=str(RESULTS / "stage3" / "specificity_survivors.json"),
    params:
        mock=MOCK,
        fixture=str(FIXTURES / "crosspan" / "sample_panel.json"),
    shell:
        r"""
        mkdir -p $(dirname {output.panel})
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture} {output.panel}
            printf '{{"mock": true, "survivors": ["design_00000"]}}\n' > {output.spec_survivors}
        else
            python workflow/scripts/crosspan.py --config {input.survivors} --out {output.panel}
        fi
        """
