rule active_learning:
    input:
        embeddings=str(RESULTS / "stage4" / "embeddings.npz"),
        fps_selection=str(RESULTS / "stage4" / "fps_selection.json"),
    output:
        predictions=str(RESULTS / "stage5" / "predictions.json"),
        next_candidates=str(RESULTS / "stage5" / "next_cycle_candidates.json"),
    params:
        mock=MOCK,
        fixture=str(FIXTURES / "active_learning" / "sample_predictions.json"),
    shell:
        r"""
        mkdir -p $(dirname {output.predictions})
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture} {output.predictions}
            printf '{{"mock": true, "next_candidates": ["design_00000"]}}\n' > {output.next_candidates}
        else
            python workflow/scripts/active_learning.py --config {input.fps_selection} --out {output.predictions}
        fi
        """
