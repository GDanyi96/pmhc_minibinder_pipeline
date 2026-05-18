rule rfdiffusion:
    input:
        target_yaml=str(RESULTS / "stage0" / "target.yaml"),
        cfg="configs/rfdiffusion_default.yaml",
    output:
        manifest=str(RESULTS / "stage1" / "manifest.json"),
        sample_pdb=str(RESULTS / "stage1" / "backbones" / "design_00000.pdb"),
    params:
        mock=MOCK,
        fixture=str(FIXTURES / "rfdiffusion" / "sample.pdb"),
    shell:
        r"""
        mkdir -p $(dirname {output.sample_pdb})
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture} {output.sample_pdb}
            printf '{{"mock": true, "num_designs": 1}}\n' > {output.manifest}
        else
            python workflow/scripts/run_rfdiffusion.py --config {input.cfg} --out {output.manifest}
        fi
        """
