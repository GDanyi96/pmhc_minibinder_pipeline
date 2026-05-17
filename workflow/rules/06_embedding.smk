rule embedding:
    input:
        spec_survivors=str(RESULTS / "stage3" / "specificity_survivors.json"),
    output:
        embeddings=str(RESULTS / "stage4" / "embeddings.npz"),
        fps_selection=str(RESULTS / "stage4" / "fps_selection.json"),
    params:
        mock=MOCK,
        fixture=str(FIXTURES / "embeddings" / "sample.npz"),
    shell:
        r"""
        mkdir -p $(dirname {output.embeddings})
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture} {output.embeddings}
            printf '{{"mock": true, "selected": ["design_00000"]}}\n' > {output.fps_selection}
        else
            python workflow/scripts/embed_designs.py --config {input.spec_survivors} --out {output.embeddings}
        fi
        """
