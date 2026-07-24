import os

configfile: "config/config.yaml"

SCRIPTS = workflow.basedir
OUT  = config["output_root"]
LOGS = config["logs_root"]

S1B = os.path.join(OUT, "Step1b_RSG_Mapping")
S2  = os.path.join(OUT, "Step2_RefCheck")
S7  = os.path.join(OUT, "Step7_Merging")
S8  = os.path.join(OUT, "Step8_Rescue")

DNA_SEQ = config["dna_seq"]
THESIS = config["thesis_dir"]
IDBASE = os.path.join(THESIS, config["idbase_subdir"])


def script(rel):
    return os.path.join(SCRIPTS, rel)


rule all:
    input:
        os.path.join(S8, "lovd_flat_final.tsv"),
        os.path.join(S8, "rescue_metrics_comparison.csv"),
        os.path.join(S8, "lovd_review_final.xlsx"),


rule classify_unresolved:
    input:
        os.path.join(S7, "unresolved_variants.tsv"),
    output:
        disp=os.path.join(S8, "unresolved_disposition.tsv"),
        summary=os.path.join(S8, "unresolved_disposition_summary.tsv"),
    conda:
        "envs/python.yaml"
    shell:
        f'IN_PATH="{{input}}" OUT_PATH="{{output.disp}}" SUMMARY="{{output.summary}}" '
        f'python "{script("Step7_Merging/classify_unresolved.py")}"'


rule generate_rescue_input:
    input:
        disp=os.path.join(S8, "unresolved_disposition.tsv"),
        lrg=os.path.join(S1B, "LRG_with_NM.csv"),
    output:
        tsv=os.path.join(S8, "rescue_mutalyzer_input.tsv"),
        skip=os.path.join(S8, "rescue_mutalyzer_skipped.tsv"),
    params:
        mane_nm=DNA_SEQ["mane_nm"],
    conda:
        "envs/python.yaml"
    shell:
        f'IN_DISPOSITION="{{input.disp}}" LRG_NM_CSV="{{input.lrg}}" '
        f'MANE_NM_DIR="{{params.mane_nm}}" '
        f'OUT_TSV="{{output.tsv}}" LOG_TSV="{{output.skip}}" '
        f'python "{script("Step7_Merging/generate_rescue_input.py")}"'


rule mutalyzer_rescue:
    input:
        os.path.join(S8, "rescue_mutalyzer_input.tsv"),
    output:
        os.path.join(S8, "rescue_mutalyzer_results.tsv"),
    params:
        cache=config["mutalyzer"]["cache_dir"],
        logs=LOGS,
        mane_nm=DNA_SEQ["mane_nm"],
        idr_nm=DNA_SEQ["idrefseq_nm"],
    conda:
        "envs/http.yaml"
    shell:
        f'INPUT_TSV="{{input}}" OUT_TSV="{{output}}" '
        f'CACHE_JSONL="{{params.cache}}/rescue_MANE.jsonl" '
        f'LOG_PATH="{{params.logs}}/mutalyzer_rescue.log" '
        f'MANE_NM_DIR="{{params.mane_nm}}" IDREFSEQ_NM_DIR="{{params.idr_nm}}" '
        f'python "{script("Step6_Mutalyzer/run_mutalyzer_NM.py")}" --source MANE'


rule vv_rescue:
    input:
        os.path.join(S8, "unresolved_disposition.tsv"),
    output:
        os.path.join(S8, "vv_rescue_results.tsv"),
    params:
        cache=config["variantvalidator"]["cache_dir"],
        logs=LOGS,
    conda:
        "envs/http.yaml"
    shell:
        f'IN_TSV="{{input}}" OUT_TSV="{{output}}" '
        f'CACHE_JSONL="{{params.cache}}/vv_rescue.jsonl" '
        f'LOG_PATH="{{params.logs}}/vv_rescue.log" '
        f'TRANSCRIPTS="mane_select" FILTER_DISPOSITIONS="RESCUABLE_AUTO" '
        f'python "{script("Step7_Merging/vv_client.py")}"'


rule resolve_unresolved:
    input:
        disp=os.path.join(S8, "unresolved_disposition.tsv"),
        mutalyzer=os.path.join(S8, "rescue_mutalyzer_results.tsv"),
        vv=os.path.join(S8, "vv_rescue_results.tsv"),
        offset=os.path.join(S2, "lrg_offset_results.csv"),
        merged=os.path.join(S7, "merged_variants.tsv"),
    output:
        resolved=os.path.join(S8, "resolved_unresolved.tsv"),
        funnel=os.path.join(S8, "resolution_funnel.tsv"),
        narrative=os.path.join(S8, "resolution_narrative.md"),
    params:
        vv_cache=config["variantvalidator"]["cache_dir"],
        logs=LOGS,
    conda:
        "envs/http.yaml"
    shell:
        f'IN_DISPOSITION="{{input.disp}}" '
        f'IN_MUTALYZER="{{input.mutalyzer}}" '
        f'IN_VV="{{input.vv}}" '
        f'OFFSET_CSV="{{input.offset}}" '
        f'MERGED_TSV="{{input.merged}}" '
        f'OUT_RESOLVED="{{output.resolved}}" '
        f'OUT_FUNNEL="{{output.funnel}}" '
        f'OUT_NARRATIVE="{{output.narrative}}" '
        f'VV_CACHE_JSONL="{{params.vv_cache}}/vv_rescue.jsonl" '
        f'LOG_PATH="{{params.logs}}/resolve_unresolved.log" '
        f'python "{script("Step7_Merging/resolve_unresolved.py")}"'


rule augment_merged:
    input:
        resolved=os.path.join(S8, "resolved_unresolved.tsv"),
        merged=os.path.join(S7, "merged_variants.tsv"),
    output:
        tsv=os.path.join(S8, "augmented_merged_variants.tsv"),
        stats=os.path.join(S8, "augment_stats.txt"),
    conda:
        "envs/python.yaml"
    shell:
        f'MERGED_TSV="{{input.merged}}" RESOLVED_TSV="{{input.resolved}}" '
        f'OUT_TSV="{{output.tsv}}" STATS_TXT="{{output.stats}}" '
        f'python "{script("Step7_Merging/merge_with_rescued.py")}"'


rule dedup_rescued:
    input:
        os.path.join(S8, "augmented_merged_variants.tsv"),
    output:
        dedup=os.path.join(S8, "dedup_merged_variants.rescued.tsv"),
        stats=os.path.join(S8, "dedup_rescued_stats.txt"),
    conda:
        "envs/python.yaml"
    shell:
        f'INPUT_TSV="{{input}}" OUTPUT_TSV="{{output.dedup}}" STATS_TXT="{{output.stats}}" '
        f'python "{script("Step7_Merging/dedup_merged_variants.py")}"'


rule build_final_lovd:
    input:
        rescued_dedup=os.path.join(S8, "dedup_merged_variants.rescued.tsv"),
        primary_dedup=os.path.join(S7, "dedup_merged_variants.tsv"),
        merged=os.path.join(S7, "merged_variants.tsv"),
        disp=os.path.join(S8, "unresolved_disposition.tsv"),
    output:
        flat=os.path.join(S8, "lovd_flat_final.tsv"),
        metrics=os.path.join(S8, "rescue_metrics_comparison.csv"),
    params:
        scripts=SCRIPTS,
    conda:
        "envs/python.yaml"
    shell:
        f'RESCUED_DEDUP_TSV="{{input.rescued_dedup}}" '
        f'PRIMARY_DEDUP_TSV="{{input.primary_dedup}}" '
        f'MERGED_TSV="{{input.merged}}" '
        f'UNRESOLVED_DISP_TSV="{{input.disp}}" '
        f'SCRIPTS_DIR="{{params.scripts}}" '
        f'OUT_FLAT="{{output.flat}}" '
        f'OUT_METRICS="{{output.metrics}}" '
        f'python "{script("Step8_Rescue/build_final_lovd.py")}"'


rule lovd_flat_final_cols:
    input:
        os.path.join(S8, "dedup_merged_variants.rescued.tsv"),
    output:
        os.path.join(S8, "lovd_flat_final_cols.tsv"),
    conda:
        "envs/python.yaml"
    shell:
        f'IN_PATH="{{input}}" OUT_PATH="{{output}}" '
        f'python "{script("Step7_Merging/build_lovd_flat.py")}"'


rule patient_metadata_final:
    input:
        flat=os.path.join(S8, "lovd_flat_final_cols.tsv"),
    output:
        os.path.join(S8, "lovd_flat_with_patients_final.tsv"),
    params:
        idbase=IDBASE,
        logs=LOGS,
    conda:
        "envs/python.yaml"
    shell:
        f'IN_PATH="{{input.flat}}" IDBASE_DIR="{{params.idbase}}" '
        f'OUT_PATH="{{output}}" LOG_PATH="{{params.logs}}/patient_metadata_missing_final.tsv" '
        f'python "{script("Step7_Merging/extract_patient_metadata.py")}"'


rule lovd_review_final:
    input:
        withpat=os.path.join(S8, "lovd_flat_with_patients_final.tsv"),
        offset=os.path.join(S2, "lrg_offset_results.csv"),
    output:
        os.path.join(S8, "lovd_review_final.xlsx"),
    conda:
        "envs/python.yaml"
    shell:
        f'IN_PATH="{{input.withpat}}" OFFSET_PATH="{{input.offset}}" OUT_PATH="{{output}}" '
        f'python "{script("Step7_Merging/build_lovd_review_xlsx.py")}"'
