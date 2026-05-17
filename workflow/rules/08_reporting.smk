rule reporting:
    input:
        metrics=str(RESULTS / "stage2" / "metrics.json"),
        panel=str(RESULTS / "stage3" / "crosspan_panel.json"),
        fps_selection=str(RESULTS / "stage4" / "fps_selection.json"),
        next_candidates=str(RESULTS / "stage5" / "next_cycle_candidates.json"),
    output:
        report=f"reports/cycle_{CYCLE}.md",
        manifest=str(RESULTS / "stage6" / "manifest.json"),
    params:
        mock=MOCK,
        fixture=str(FIXTURES / "report" / "sample_manifest.json"),
        cycle=CYCLE,
    shell:
        r"""
        mkdir -p $(dirname {output.manifest}) reports
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture} {output.manifest}
            printf '# Cycle %s — mock report\n\nGenerated in mock mode. See `INDEX.md` and `specs/stage6_reporting.md`.\n' "{params.cycle}" > {output.report}
        else
            python workflow/scripts/render_report.py --config {input.metrics} --out {output.report}
        fi
        """
