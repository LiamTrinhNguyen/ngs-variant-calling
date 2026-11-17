#!/usr/bin/env nextflow
/*
 * =============================================================================
 * Nextflow pipeline – NGS variant calling for many samples
 * =============================================================================
 * WHAT: the same steps as run_interactive_pipeline.py and workflow/Snakefile,
 *       written as Nextflow DSL2 processes. Reads its samples and filters from
 *       the SAME config/config.yaml as Snakemake, and writes the same numbered
 *       folders (00_raw … 07_dashboard) plus pipeline_info/ (Nextflow's run
 *       report, timeline and trace).
 * WHY:  Nextflow is the other widely used workflow engine (nf-core). It runs the
 *       same code unchanged on a laptop, an HPC scheduler (SLURM, PBS) or the
 *       cloud (AWS Batch, Google Batch) by changing only the `executor` in
 *       nextflow.config.
 *
 * Each process calls a function from run_interactive_pipeline.py through a
 * short Python snippet.
 * WHY: the commands and filters are defined in one place, so all three ways of
 *      running give the same results. If a function changes, update the
 *      matching process here AND the rule in workflow/Snakefile.
 *
 * Usage (from ngs-variant-calling/, with ngs-env active):
 *   nextflow run workflow/main.nf                               # config/config.yaml
 *   nextflow run workflow/main.nf --config my.yaml --outdir /path/to/out
 *   nextflow run workflow/main.nf -resume                        # reuse finished steps
 * =============================================================================
 */

nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// Python header used by every process: import the pipeline module.
// WHY projectDir/..: main.nf is in workflow/, the module one folder up.
// -----------------------------------------------------------------------------
def py(String body) {
    """
    python3 - <<'PY'
    import sys, shutil
    from pathlib import Path
    sys.path.insert(0, "${projectDir.parent}")
    import run_interactive_pipeline as ngs
    ${body.stripIndent().trim().readLines().join('\n    ')}
    PY
    """.stripIndent()
}

// The reference FASTA among the staged reference files (the rest are its indexes)
def fastaOf(files) { files.find { it.name.endsWith('.fasta') } }

// =============================================================================
// Processes – one per step. publishDir copies each step's files into the same
// numbered folder the other two drivers use.
// WHY mode 'copy': the work/ folder can be deleted later without losing results.
// =============================================================================

// STEP 2 · reference genome (NCBI), shared by all samples
process DOWNLOAD_REFERENCE {
    tag "${ref}"
    publishDir "${params.outdir}/01_reference", mode: 'copy'

    input:
    val ref

    output:
    path "${ref}.fasta"

    script:
    py("""
    ngs.download_reference("unknown", Path("."), force_accession="${ref}")
    """)
}

// STEP 3 · indexing (bwa index, samtools faidx/dict), built once for all samples
process INDEX_REFERENCE {
    tag "${fasta.baseName}"
    publishDir "${params.outdir}/01_reference", mode: 'copy', pattern: "*.{amb,ann,bwt,pac,sa,fai,dict}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/*"

    input:
    path fasta

    output:
    // WHY includeInputs: later steps need the FASTA and its index files together
    path "${fasta.baseName}*", includeInputs: true, emit: ref
    path "logs/*", emit: logs

    script:
    py("""
    Path("logs").mkdir(exist_ok=True)
    ngs.index_reference(Path("${fasta}"), Path("logs"))
    """)
}

// STEP 4 · data input – SRA samples only (local r1:/r2: samples skip this)
process DOWNLOAD_READS {
    tag "${sample}"
    publishDir "${params.outdir}/00_raw", mode: 'copy'

    input:
    tuple val(sample), val(acc)

    output:
    tuple val(sample), path("${sample}_1.fastq.gz"), path("${sample}_2.fastq.gz")

    script:
    // WHY a private download folder: ENA names files after the accession, so
    //     two samples with the same SRR must not share one folder
    py("""
    tmp = Path("download"); tmp.mkdir()
    _, fastqs = ngs.get_sra_metadata_and_fastq("${acc}", tmp)
    if len(fastqs) != 2:
        sys.exit("${acc} is not paired-end; this workflow needs paired reads")
    shutil.move(str(fastqs[0]), "${sample}_1.fastq.gz")
    shutil.move(str(fastqs[1]), "${sample}_2.fastq.gz")
    shutil.rmtree(tmp)
    """)
}

// STEP 5 · QC & trimming (fastp)
process FASTP {
    tag "${sample}"
    publishDir "${params.outdir}/02_qc", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_R1.trimmed.fastq.gz"), path("${sample}_R2.trimmed.fastq.gz"), emit: reads
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.{html,log}"

    script:
    py("""
    ngs.run_fastp(Path("${r1}"), Path("${r2}"), Path("."), "${sample}", ${task.cpus})
    """)
}

// STEP 6 · alignment (BWA-MEM + samtools sort/index)
process ALIGN {
    tag "${sample}"
    publishDir "${params.outdir}/03_aligned", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)
    path ref_files

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")

    script:
    py("""
    ngs.align_and_sort(Path("${fastaOf(ref_files)}"), Path("${r1}"), Path("${r2}"),
                       Path("."), "${sample}", ${task.cpus})
    """)
}

// STEP 7 · variant calling (bcftools mpileup | call --ploidy 1 | norm)
process CALL_VARIANTS {
    tag "${sample}"
    publishDir "${params.outdir}/04_variants", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_files

    output:
    tuple val(sample), path("${sample}.raw.vcf.gz"), path("${sample}.raw.vcf.gz.csi"), emit: vcf
    path "${sample}.calling_settings.json", emit: settings   // read by the dashboard
    path "${sample}.bcftools.log"

    script:
    py("""
    ngs.call_variants(Path("${fastaOf(ref_files)}"), Path("${bam}"), Path("."), "${sample}")
    """)
}

// STEP 8 · annotation (SnpEff) – only used when SNPEFF_DB has the reference
// NOTE: like Snakemake, a failure STOPS the run (the interactive script continues).
process ANNOTATE {
    tag "${sample}"
    publishDir "${params.outdir}/05_annotated", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(csi)
    val genome

    output:
    tuple val(sample), path("${sample}.ann.vcf.gz"), path("${sample}.ann.vcf.gz.tbi")

    script:
    py("""
    ngs.annotate_variants(Path("${vcf}"), "${genome}", Path("."), "${sample}")
    """)
}

// STEP 9 · flag & report (bcftools filter + Polars)
process REPORT {
    tag "${sample}"
    publishDir "${params.outdir}/06_reports", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(idx)
    val opts      // filters: min_qual, min_dp, min_af

    output:
    tuple val(sample), path("${sample}.all_calls.tsv"), emit: calls
    path "${sample}.high_confidence_variants.tsv", emit: report
    path "${sample}.flagged.vcf.gz*"

    script:
    py("""
    ngs.filter_and_report(Path("${vcf}"), Path("."), "${sample}",
                          min_qual=${opts.min_qual}, min_dp=${opts.min_dp}, min_af=${opts.min_af})
    """)
}

// STEP 9b · statistics in R (R/variant_stats.R, Rscript from the ngs-r env)
process R_STATS {
    tag "${sample}"
    publishDir "${params.outdir}/06_reports", mode: 'copy'

    input:
    tuple val(sample), path(calls)
    path ref_files
    val opts      // hotspot_window

    output:
    path "${sample}.r_{stats,windows}.tsv", emit: tables
    path "${sample}.r_stats.log"

    script:
    py("""
    ngs.run_r_stats(Path("${calls}"), Path("${fastaOf(ref_files)}"), Path("."),
                    "${sample}", ${opts.hotspot_window})
    """)
}

// STEP 10 · one dashboard for all samples
// WHAT: stages every sample's files into the numbered layout that
//       dashboard.build_dashboard() reads, then builds 07_dashboard/.
// WHY stageAs: Nextflow gives each process its own folder; stageAs recreates
//     01_reference/, 02_qc/ … inside it so the dashboard code needs no changes.
process DASHBOARD {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ref_files, stageAs: '01_reference/*'
    path qc_json,   stageAs: '02_qc/*'
    path bams,      stageAs: '03_aligned/*'
    path settings,  stageAs: '04_variants/*'
    path reports,   stageAs: '06_reports/*'
    val samples
    val opts

    output:
    path "07_dashboard/*"

    script:
    def names = samples.collect { "\"${it}\"" }.join(', ')
    py("""
    import dashboard
    dashboard.build_dashboard(Path("."), min_qual=${opts.min_qual}, min_dp=${opts.min_dp},
                              min_af=${opts.min_af}, samples=[${names}])
    """)
}

// =============================================================================
// Workflow – connects the processes with channels
// =============================================================================
workflow {
    // ---- settings: config/config.yaml (shared with Snakemake) ----
    // WHY the same file: one list of samples and filters for both engines.
    //     Command-line --reference / --min_qual … still override it (see nextflow.config).
    // WHY a map (opts) instead of changing params: params are read-only once
    //     the workflow starts, so the merged settings are passed to the processes.
    def cfg = new org.yaml.snakeyaml.Yaml().load(file(params.config).text)
    def opts = [
        reference:      params.reference      ?: cfg.reference,
        min_qual:       params.min_qual       ?: (cfg.min_qual ?: 30),
        min_dp:         params.min_dp         ?: (cfg.min_dp ?: 15),
        min_af:         params.min_af         ?: (cfg.min_af ?: 0.8),
        hotspot_window: params.hotspot_window ?: (cfg.hotspot_window ?: 1000),
    ]
    def samples = cfg.samples as Map
    log.info "NGS variant calling (Nextflow) · ${samples.size()} sample(s) · reference ${opts.reference} · " +
             "filters QUAL>=${opts.min_qual} DP>=${opts.min_dp} AF>=${opts.min_af} · out ${params.outdir}"

    // Ask the Python module whether SnpEff has a database for this reference.
    // WHY ask Python: SNPEFF_DB is defined once, in run_interactive_pipeline.py.
    def probe = ["python3", "-c",
        "import sys; sys.path.insert(0, '${projectDir.parent}'); " +
        "import run_interactive_pipeline as m; print(m.SNPEFF_DB.get('${opts.reference}') or '')"].execute()
    probe.waitFor()
    def snpeff_genome = probe.text.trim()
    if (!snpeff_genome)
        log.warn "Reference ${opts.reference} has no entry in SNPEFF_DB – annotation will be skipped."

    // ---- reference: download + index once ----
    ref_files = INDEX_REFERENCE(DOWNLOAD_REFERENCE(opts.reference)).ref.collect()

    // ---- reads: SRA samples are downloaded, local samples used in place ----
    sra_ch   = channel.fromList(samples.findAll { k, v -> v.sra }.collect { k, v -> tuple(k, v.sra) })
    local_ch = channel.fromList(samples.findAll { k, v -> v.r1 && v.r2 }
                                       .collect { k, v -> tuple(k, file(v.r1), file(v.r2)) })
    reads_ch = DOWNLOAD_READS(sra_ch).mix(local_ch)

    // ---- per-sample steps ----
    FASTP(reads_ch)
    bam_ch = ALIGN(FASTP.out.reads, ref_files)
    CALL_VARIANTS(bam_ch, ref_files)
    vcf_ch = snpeff_genome ? ANNOTATE(CALL_VARIANTS.out.vcf, snpeff_genome) : CALL_VARIANTS.out.vcf
    REPORT(vcf_ch, opts)
    R_STATS(REPORT.out.calls, ref_files, opts)

    // ---- one dashboard for all samples ----
    // .collect() waits until every sample has finished the step
    DASHBOARD(
        ref_files,
        FASTP.out.json.collect(),
        bam_ch.map { s, bam, bai -> [bam, bai] }.flatten().collect(),
        CALL_VARIANTS.out.settings.collect(),
        REPORT.out.report.mix(REPORT.out.calls.map { s, f -> f }, R_STATS.out.tables.flatten()).collect(),
        samples.keySet() as List,
        opts,
    )
}
