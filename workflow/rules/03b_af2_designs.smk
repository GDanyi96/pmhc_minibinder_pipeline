# Stage 2 designs marker rule -- pure DAG anchor downstream of 02b.
#
# stage2_summary.json is produced by `stage2_designs` (02b); this rule
# materializes a `.done` marker so future stages (Stage 3 cross-pan and
# beyond) can declare a single, stable upstream dependency without
# coupling to the summary file's name or schema.

rule stage2_designs_done:
    input:
        summary=str(RESULTS / "stage2" / "stage2_summary.json"),
    output:
        marker=str(RESULTS / "stage2" / ".done"),
    shell:
        r"""
        mkdir -p $(dirname {output.marker})
        touch {output.marker}
        """
