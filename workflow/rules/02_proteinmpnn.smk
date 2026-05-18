PMPNN_TARGETS = ("3hpj_clean", "2bnr_clean")


rule proteinmpnn:
    input:
        primary_pdb="data/targets/3hpj_clean.pdb",
        control_pdb="data/targets/2bnr_clean.pdb",
        chain_jsonl="configs/proteinmpnn_chains.json",
        # DAG ordering anchor: Stage 1 must finish before Stage 2 sequence
        # design starts (cycle 2 onward). Not yet consumed by the
        # subprocess call -- splicing is run_stage2_designs.py's job (next
        # session).
        stage1_summary=str(RESULTS / "stage1" / "rfdiffusion" / "stage1_summary.json"),
    output:
        primary_fasta=str(RESULTS / "stage2" / "proteinmpnn" / "3hpj_clean.fasta"),
        control_fasta=str(RESULTS / "stage2" / "proteinmpnn" / "2bnr_clean.fasta"),
    params:
        mock=MOCK,
        fixture=str(FIXTURES / "proteinmpnn" / "sample.fasta"),
        out_dir=str(RESULTS / "stage2" / "proteinmpnn"),
    shell:
        r"""
        mkdir -p {params.out_dir}
        if [ "{params.mock}" = "True" ]; then
            cp {params.fixture} {output.primary_fasta}
            cp {params.fixture} {output.control_fasta}
        else
            python workflow/scripts/run_proteinmpnn.py \
                --pdb_path {input.primary_pdb} \
                --chain_id_jsonl {input.chain_jsonl} \
                --fixed_chains A B C --designed_chain D \
                --out_folder {params.out_dir}
            python workflow/scripts/run_proteinmpnn.py \
                --pdb_path {input.control_pdb} \
                --chain_id_jsonl {input.chain_jsonl} \
                --fixed_chains A B C --designed_chain D \
                --out_folder {params.out_dir}
        fi
        """
