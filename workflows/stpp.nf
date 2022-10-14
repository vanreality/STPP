/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate input parameters
WorkflowSTPP.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.bwa,
    params.dbnsfp,
    params.dbnsfp_tbi,
    params.dbsnp,
    params.dbsnp_tbi,
    params.dict,
    params.fasta,
    params.fasta_fai,
    params.germline_resource,
    params.germline_resource_tbi,
    params.input,
    params.intervals,
    params.known_snps,
    params.known_snps_tbi,
    params.known_indels,
    params.known_indels_tbi,
    params.pon,
    params.pon_tbi,
    params.snpeff_cache,
    params.spliceai_indel,
    params.spliceai_indel_tbi,
    params.spliceai_snv,
    params.spliceai_snv_tbi,
    params.vep_cache
]

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

// Set input, can either be from --input or from automatic retrieval in WorkflowSarek.groovy
ch_input_sample = extract_csv(file(params.input, checkIfExists: true))

// Fails when no intervals file provided
if (!params.intervals){
    log.error "Target file specified with `--intervals` must be provided"
    exit 1
}

// Fails when wrongfull extension for intervals file
if (params.intervals && !params.intervals.endsWith("bed")){
    log.error "Target file specified with `--intervals` must be in BED format for targeted data"
    exit 1
}

// Warns when missing files or params for mutect2
if(params.tools && params.tools.split(',').contains('mutect2')){
    if(!params.pon){
        log.warn("No Panel-of-normal was specified for Mutect2.\nIt is highly recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2\nFor more information on how to create one: https://gatk.broadinstitute.org/hc/en-us/articles/5358921041947-CreateSomaticPanelOfNormals-BETA-")
    }
    if(!params.germline_resource){
        log.warn("If Mutect2 is specified without a germline resource, no filtering will be done.\nIt is recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2")
    }
    // if(params.pon && params.pon.contains("/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz")){
    //     log.warn("The default Panel-of-Normals provided by GATK is used for Mutect2.\nIt is highly recommended to generate one from normal samples that are technical similar to the tumor ones.\nFor more information: https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-")
    // }
}

// Fails when missing resources for baserecalibrator
if(!params.dbsnp && !params.known_indels){
    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate'] && (!params.skip_tools || (params.skip_tools && !params.skip_tools.split(',').contains('baserecalibrator')))){
        log.error "Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command."
        exit 1
    }
}

// Fails when missing tools for variant_calling or annotate
if ((params.step == 'variant_calling' || params.step == 'annotate') && !params.tools) {
    log.error "Please specify at least one tool when using `--step ${params.step}`."
    exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
dbsnp              = params.dbsnp              ? Channel.fromPath(params.dbsnp).collect()                    : Channel.value([])
known_snps         = params.known_snps         ? Channel.fromPath(params.known_snps).collect()               : Channel.value([])
fasta              = params.fasta              ? Channel.fromPath(params.fasta).collect()                    : Channel.empty()
fasta_fai          = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()                : Channel.empty()
germline_resource  = params.germline_resource  ? Channel.fromPath(params.germline_resource).collect()        : Channel.value([]) //Mutec2 does not require a germline resource, so set to optional input
known_indels       = params.known_indels       ? Channel.fromPath(params.known_indels).collect()             : Channel.value([])
known_snps         = params.known_snps         ? Channel.fromPath(params.known_snps).collect()               : Channel.value([])
pon                = params.pon                ? Channel.fromPath(params.pon).collect()                      : Channel.value([]) //PON is optional for Mutect2 (but highly recommended)

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
snpeff_db          = params.snpeff_db          ?: Channel.empty()
vep_cache_version  = params.vep_cache_version  ?: Channel.empty()
vep_genome         = params.vep_genome         ?: Channel.empty()
vep_species        = params.vep_species        ?: Channel.empty()

// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
snpeff_cache       = params.snpeff_cache       ? Channel.fromPath(params.snpeff_cache).collect()             : []
vep_cache          = params.vep_cache          ? Channel.fromPath(params.vep_cache).collect()                : []

vep_extra_files = []

if (params.dbnsfp && params.dbnsfp_tbi) {
    vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
    vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Create samplesheets to restart from different steps
include { MAPPING_CSV                                          } from '../subworkflows/local/mapping_csv'
include { MARKDUPLICATES_CSV                                   } from '../subworkflows/local/markduplicates_csv'
include { PREPARE_RECALIBRATION_CSV                            } from '../subworkflows/local/prepare_recalibration_csv'
include { RECALIBRATE_CSV                                      } from '../subworkflows/local/recalibrate_csv'
include { VARIANTCALLING_CSV                                   } from '../subworkflows/local/variantcalling_csv'

// Build indices if needed
include { PREPARE_GENOME                                       } from '../subworkflows/local/prepare_genome'

// Build intervals if needed
include { PREPARE_INTERVALS                                    } from '../subworkflows/local/prepare_intervals'

// Run FASTQC
include { RUN_FASTQC                                           } from '../subworkflows/nf-core/run_fastqc'

// TRIM/SPLIT FASTQ Files
include { FASTP                                                } from '../modules/fastp/main'

// Map input reads to reference genome
include { GATK4_MAPPING                                        } from '../subworkflows/nf-core/gatk4/mapping/main'

// Merge and index BAM files (optional)
include { MERGE_INDEX_BAM                                      } from '../subworkflows/nf-core/merge_index_bam'

// Convert bam to cram and cram to bam
include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM               } from '../modules/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM_RECAL         } from '../modules/samtools/convert/main'

include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM               } from '../modules/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM_VARIANTCALLING} from '../modules/samtools/convert/main'

// Mark Duplicates (+QC)
include { MARKDUPLICATES                                       } from '../subworkflows/nf-core/gatk4/markduplicates/main'

// Convert to CRAM (+QC)
include { BAM_TO_CRAM                                          } from '../subworkflows/nf-core/bam_to_cram'

// Create recalibration tables
include { PREPARE_RECALIBRATION                                } from '../subworkflows/nf-core/gatk4/prepare_recalibration/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { RECALIBRATE                                          } from '../subworkflows/nf-core/gatk4/recalibrate/main'

// QC on CRAM
include { CRAM_QC                                              } from '../subworkflows/nf-core/cram_qc'

// Variant calling on tumor/normal pair
include { PAIR_VARIANT_CALLING                                 } from '../subworkflows/local/pair_variant_calling'

// Variant results QC
include { VCF_QC                                               } from '../subworkflows/nf-core/vcf_qc'

// Annotation
include { ANNOTATE                                             } from '../subworkflows/local/annotate'










/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STPP{

    // To gather all QC reports for MultiQC
    ch_reports  = Channel.empty()
    // To gather used softwares versions for MultiQC
    ch_versions = Channel.empty()


    // Build indices if needed
    PREPARE_GENOME(
        dbsnp,
        fasta,
        fasta_fai,
        germline_resource,
        known_indels,
        known_snps,
        pon)

    // Gather built indices or get them from the params
    bwa                    = params.fasta                   ? params.bwa                        ? Channel.fromPath(params.bwa).collect()                   : PREPARE_GENOME.out.bwa                   : []
    dict                   = params.fasta                   ? params.dict                       ? Channel.fromPath(params.dict).collect()                  : PREPARE_GENOME.out.dict                  : []
    fasta_fai              = params.fasta                   ? params.fasta_fai                  ? Channel.fromPath(params.fasta_fai).collect()             : PREPARE_GENOME.out.fasta_fai             : []
    dbsnp_tbi              = params.dbsnp                   ? params.dbsnp_tbi                  ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_GENOME.out.dbsnp_tbi             : Channel.value([])
    germline_resource_tbi  = params.germline_resource       ? params.germline_resource_tbi      ? Channel.fromPath(params.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : []
    known_snps_tbi         = params.known_snps              ? params.known_snps_tbi             ? Channel.fromPath(params.known_snps_tbi).collect()        : PREPARE_GENOME.out.known_snps_tbi        : Channel.value([])
    known_indels_tbi       = params.known_indels            ? params.known_indels_tbi           ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_GENOME.out.known_indels_tbi      : Channel.value([])
    pon_tbi                = params.pon                     ? params.pon_tbi                    ? Channel.fromPath(params.pon_tbi).collect()               : PREPARE_GENOME.out.pon_tbi               : []
    msisensorpro_scan      = PREPARE_GENOME.out.msisensorpro_scan

    // Gather index for mapping given the chosen aligner
    ch_map_index = bwa

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels     = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    known_sites_snps     = dbsnp.concat(known_snps).collect()
    known_sites_snps_tbi = dbsnp_tbi.concat(known_snps_tbi).collect()

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    intervals_bed_combined      = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined  // [interval.bed] all intervals in one file
    intervals_for_preprocessing = []      // For QC during preprocessing, we don't need any intervals

    intervals                   = PREPARE_INTERVALS.out.intervals_bed        // [interval, num_intervals] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi        = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [interval_bed, tbi, num_intervals] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather

    // Gather used softwares versions
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // STEP 1: MAPPING
    if (params.step == 'mapping') {

        // BEFORE MAPPING: QC & TRIM
        // `--skip_tools fastqc` to skip fastqc
        // trim only with `--trim_fastq`
        // additional options to be set up

        // QC
        if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
            RUN_FASTQC(ch_input_sample)

            ch_reports  = ch_reports.mix(RUN_FASTQC.out.fastqc_zip.collect{meta, logs -> logs})
            ch_versions = ch_versions.mix(RUN_FASTQC.out.versions)
        }

        ch_reads_fastp = ch_input_sample

        // Trimming and/or splitting
        if (params.trim_fastq) {

            save_trimmed_fail = false
            save_merged = false
            FASTP(ch_reads_fastp, save_trimmed_fail, save_merged)

            ch_reports = ch_reports.mix(
                                    FASTP.out.json.collect{meta, json -> json},
                                    FASTP.out.html.collect{meta, html -> html}
                                    )

            if(params.split_fastq){
                ch_reads_to_map = FASTP.out.reads.map{ key, reads ->

                        read_files = reads.sort{ a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                        [[
                            data_type:key.data_type,
                            id:key.id,
                            numLanes:key.numLanes,
                            patient: key.patient,
                            read_group:key.read_group,
                            sample:key.sample,
                            sex:key.sex,
                            size:read_files.size(),
                            status:key.status,
                        ],
                        read_files]
                    }.transpose()
            }else{
                ch_reads_to_map = FASTP.out.reads
            }

            ch_versions = ch_versions.mix(FASTP.out.versions)
        } else {
            ch_reads_to_map = ch_reads_fastp
        }

        // MAPPING READS TO REFERENCE GENOME
        // reads will be sorted
        ch_reads_to_map = ch_reads_to_map.map{ meta, reads ->
            // update ID when no multiple lanes or splitted fastqs
            new_id = meta.size * meta.numLanes == 1 ? meta.sample : meta.id

            [[
                data_type:  meta.data_type,
                id:         new_id,
                numLanes:   meta.numLanes,
                patient:    meta.patient,
                read_group: meta.read_group,
                sample:     meta.sample,
                sex:        meta.sex,
                size:       meta.size,
                status:     meta.status,
                ],
            reads]
        }

        sort_bam = true
        GATK4_MAPPING(ch_reads_to_map, ch_map_index, sort_bam)

        // Grouping the bams from the same samples not to stall the workflow
        ch_bam_mapped = GATK4_MAPPING.out.bam.map{ meta, bam ->
            numLanes = meta.numLanes ?: 1
            size     = meta.size     ?: 1

            // update ID to be based on the sample name
            // update data_type
            // remove no longer necessary fields:
            //   read_group: Now in the BAM header
            //     numLanes: Was only needed for mapping
            //         size: Was only needed for mapping
            new_meta = [
                        id:meta.sample,
                        data_type:"bam",
                        patient:meta.patient,
                        sample:meta.sample,
                        sex:meta.sex,
                        status:meta.status,
                    ]

            // Use groupKey to make sure that the correct group can advance as soon as it is complete
            // and not stall the workflow until all reads from all channels are mapped
            [ groupKey(new_meta, numLanes * size), bam]
        }.groupTuple()

        // gatk4 markduplicates can handle multiple bams as input, so no need to merge/index here
        // Except if and only if skipping markduplicates or saving mapped bams
        if (params.save_bam_mapped || (params.skip_tools && params.skip_tools.split(',').contains('markduplicates'))) {

            // bams are merged (when multiple lanes from the same sample), indexed and then converted to cram
            MERGE_INDEX_BAM(ch_bam_mapped)

            // Create CSV to restart from this step
            MAPPING_CSV(MERGE_INDEX_BAM.out.bam_bai)

            // Gather used softwares versions
            ch_versions = ch_versions.mix(MERGE_INDEX_BAM.out.versions)
        }

        // Gather used softwares versions
        ch_versions = ch_versions.mix(GATK4_MAPPING.out.versions)
    }

    // STEP 2: markduplicates (+QC) + convert to CRAM
    if (params.step in ['mapping', 'markduplicates']) {

        // ch_bam_for_markduplicates will countain bam mapped with GATK4_MAPPING when step is mapping
        // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
        // ch_bam_for_markduplicates = params.step == 'mapping'? ch_bam_mapped : ch_input_sample.map{ meta, input, index -> [meta, input] }

        ch_bam_for_markduplicates   = Channel.empty()
        ch_input_cram_indexed       = Channel.empty()
        ch_cram_no_markduplicates   = Channel.empty()
        ch_cram_markduplicates      = Channel.empty()

        if (params.step == 'mapping') ch_bam_for_markduplicates = ch_bam_mapped
        else {
            ch_input_sample.branch{
                bam:  it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }.set{ch_convert}

            ch_bam_for_markduplicates = ch_convert.bam.map{ meta, bam, bai -> [meta, bam]}

            //In case Markduplicates is run convert CRAM files to BAM, because the tool only runs on BAM files.
            if (!(params.skip_tools && params.skip_tools.split(',').contains('markduplicates'))){

                SAMTOOLS_CRAMTOBAM(ch_convert.cram, fasta, fasta_fai)
                ch_versions = ch_versions.mix(SAMTOOLS_CRAMTOBAM.out.versions)

                ch_bam_for_markduplicates = ch_bam_for_markduplicates.mix(SAMTOOLS_CRAMTOBAM.out.alignment_index.map{ meta, bam, bai -> [meta, bam]})
            } else {
                ch_input_cram_indexed     = ch_convert.cram
            }
        }

        if (params.skip_tools && params.skip_tools.split(',').contains('markduplicates')) {

            // ch_bam_indexed will countain bam mapped with GATK4_MAPPING when step is mapping
            // which are then merged and indexed
            // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
            ch_bam_indexed = params.step == 'mapping' ? MERGE_INDEX_BAM.out.bam_bai : ch_convert.bam

            BAM_TO_CRAM(
                ch_bam_indexed,
                ch_input_cram_indexed,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            ch_cram_no_markduplicates = BAM_TO_CRAM.out.cram_converted

            // Gather QC reports
            ch_reports  = ch_reports.mix(BAM_TO_CRAM.out.qc.collect{meta, report -> report})

            // Gather used softwares versions
            ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions)
        } else {
            MARKDUPLICATES(
                ch_bam_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            ch_cram_markduplicates = MARKDUPLICATES.out.cram

            // Gather QC reports
            ch_reports  = ch_reports.mix(MARKDUPLICATES.out.qc.collect{meta, report -> report})

            // Gather used softwares versions
            ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions)
        }

        // ch_md_cram_for_restart contains either:
        // - crams from markduplicates
        // - crams converted from bam mapped when skipping markduplicates
        ch_md_cram_for_restart = Channel.empty().mix(
            ch_cram_markduplicates,
            ch_cram_no_markduplicates).map{ meta, cram, crai ->
                        //Make sure correct data types are carried through
                        [[
                            data_type:  "cram",
                            id:         meta.id,
                            patient:    meta.patient,
                            sample:     meta.sample,
                            sex:        meta.sex,
                            status:     meta.status
                            ],
                        cram, crai]
                    }

        // CSV should be written for the file actually out, either CRAM or BAM
        // Create CSV to restart from this step
        if (!(params.skip_tools && params.skip_tools.split(',').contains('markduplicates'))) MARKDUPLICATES_CSV(ch_md_cram_for_restart)
    }

    // STEP 3: prepare BQSR
    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration']) {

        // Run if starting from step "prepare_recalibration" to get input
        if(params.step == 'prepare_recalibration'){

            //Support if starting from BAM or CRAM files
            ch_input_sample.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }.set{ch_convert}

            //BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            SAMTOOLS_BAMTOCRAM(ch_convert.bam, fasta, fasta_fai)
            ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM.out.versions)

            ch_cram_for_prepare_recalibration = Channel.empty().mix(SAMTOOLS_BAMTOCRAM.out.alignment_index, ch_convert.cram)

            ch_md_cram_for_restart = SAMTOOLS_BAMTOCRAM.out.alignment_index

        } else {

            // ch_cram_for_prepare_recalibration contains either:
            // - crams from markduplicates
            // - crams converted from bam mapped when skipping markduplicates
            // - input cram files, when start from step markduplicates
            //ch_md_cram_for_restart.view() //contains md.cram.crai
            ch_cram_for_prepare_recalibration = Channel.empty().mix(ch_md_cram_for_restart, ch_input_cram_indexed)
        }

        // Create recalibration tables
        if (!(params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'))) {
            ch_table_bqsr = Channel.empty()

            PREPARE_RECALIBRATION(
                ch_cram_for_prepare_recalibration,
                dict,
                fasta,
                fasta_fai,
                intervals,
                known_sites_indels,
                known_sites_indels_tbi)

            ch_table_bqsr = PREPARE_RECALIBRATION.out.table_bqsr

            // Gather used softwares versions
            ch_versions = ch_versions.mix(PREPARE_RECALIBRATION.out.versions)

            ch_reports  = ch_reports.mix(ch_table_bqsr.collect{ meta, table -> table})

            ch_cram_applybqsr = ch_cram_for_prepare_recalibration.join(ch_table_bqsr)

            // Create CSV to restart from this step
            PREPARE_RECALIBRATION_CSV(ch_md_cram_for_restart.join(ch_table_bqsr), params.skip_tools)
        }
    }

    // STEP 4: BQSR
    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate']) {

        // Run if starting from step "prepare_recalibration"
        if(params.step == 'recalibrate'){

            //Support if starting from BAM or CRAM files
            ch_input_sample.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }.set{ch_convert}

            //If BAM file, split up table and mapped file to convert BAM to CRAM
            ch_bam_table = ch_convert.bam.map{ meta, bam, bai, table -> [meta, table]}
            ch_bam_bam   = ch_convert.bam.map{ meta, bam, bai, table -> [meta, bam, bai]}

            //BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            SAMTOOLS_BAMTOCRAM(ch_bam_bam, fasta, fasta_fai)
            ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM.out.versions)

            ch_cram_applybqsr = Channel.empty().mix(
                                    SAMTOOLS_BAMTOCRAM.out.alignment_index.join(ch_bam_table),
                                    ch_convert.cram) // Join together converted cram with input tables
        }

        if (!(params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'))) {
            ch_cram_variant_calling = Channel.empty()

            RECALIBRATE(
                ch_cram_applybqsr,
                dict,
                fasta,
                fasta_fai,
                intervals)

            ch_cram_variant_calling = RECALIBRATE.out.cram

            // Gather used softwares versions
            ch_versions = ch_versions.mix(RECALIBRATE.out.versions)

            CRAM_QC(
                ch_cram_variant_calling,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            // Gather QC reports
            ch_reports  = ch_reports.mix(CRAM_QC.out.qc.collect{meta, report -> report})

            // Gather used softwares versions
            ch_versions = ch_versions.mix(CRAM_QC.out.versions)

            //If params.save_output_as_bam, then convert CRAM files to BAM
            SAMTOOLS_CRAMTOBAM_RECAL(ch_cram_variant_calling, fasta, fasta_fai)
            ch_versions = ch_versions.mix(SAMTOOLS_CRAMTOBAM_RECAL.out.versions)

            // CSV should be written for the file actually out out, either CRAM or BAM
            csv_recalibration = Channel.empty()
            csv_recalibration = params.save_output_as_bam ?  SAMTOOLS_CRAMTOBAM_RECAL.out.alignment_index : ch_cram_variant_calling

            // Create CSV to restart from this step
            RECALIBRATE_CSV(csv_recalibration)

        } else if (params.step == 'recalibrate'){
            // ch_cram_variant_calling contains either:
            // - input bams converted to crams, if started from step recal + skip BQSR
            // - input crams if started from step recal + skip BQSR
            ch_cram_variant_calling = Channel.empty().mix(SAMTOOLS_BAMTOCRAM.out.alignment_index,
                                                        ch_convert.cram.map{ meta, cram, crai, table -> [meta, cram, crai]})
        } else {
            // ch_cram_variant_calling contains either:
            // - crams from markduplicates = ch_cram_for_prepare_recalibration if skip BQSR but not started from step recalibration
            ch_cram_variant_calling = Channel.empty().mix(ch_cram_for_prepare_recalibration)
        }
    }

    // STEP 5: Variant calling and annotating
    if (params.step == 'variant_calling') {

        ch_input_sample.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }.set{ch_convert}

        //BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
        SAMTOOLS_BAMTOCRAM_VARIANTCALLING(ch_convert.bam, fasta, fasta_fai)
        ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM_VARIANTCALLING.out.versions)

        ch_cram_variant_calling = Channel.empty().mix(SAMTOOLS_BAMTOCRAM_VARIANTCALLING.out.alignment_index, ch_convert.cram)
    }

    if (params.tools) {

        if (params.step != 'annotate') {

            //
            // Logic to separate germline samples, tumor samples with no matched normal, and combine tumor-normal pairs
            //
            ch_cram_variant_calling.branch{
                normal: it[0].status == 0
                tumor:  it[0].status == 1
            }.set{ch_cram_variant_calling_status}

            // All Germline samples
            ch_cram_variant_calling_normal_to_cross = ch_cram_variant_calling_status.normal.map{ meta, cram, crai -> [meta.patient, meta, cram, crai] }

            // All tumor samples
            ch_cram_variant_calling_pair_to_cross = ch_cram_variant_calling_status.tumor.map{ meta, cram, crai -> [meta.patient, meta, cram, crai] }

            // Tumor only samples
            // 1. Group together all tumor samples by patient ID [patient1, [meta1, meta2], [cram1,crai1, cram2, crai2]]

            // Downside: this only works by waiting for all tumor samples to finish preprocessing, since no group size is provided
            ch_cram_variant_calling_tumor_grouped = ch_cram_variant_calling_pair_to_cross.groupTuple()

            // 2. Join with normal samples, in each channel there is one key per patient now. Patients without matched normal end up with: [patient1, [meta1, meta2], [cram1,crai1, cram2, crai2], null]
            ch_cram_variant_calling_tumor_joined = ch_cram_variant_calling_tumor_grouped.join(ch_cram_variant_calling_normal_to_cross, remainder: true)

            // 3. Filter out entries with last entry null
            ch_cram_variant_calling_tumor_filtered = ch_cram_variant_calling_tumor_joined.filter{ it ->  !(it.last()) }

            // 4. Transpose [patient1, [meta1, meta2], [cram1,crai1, cram2, crai2]] back to [patient1, meta1, [cram1,crai1], null] [patient1, meta2, [cram2,crai2], null]
            // and remove patient ID field & null value for further processing [meta1, [cram1,crai1]] [meta2, [cram2,crai2]]
            ch_cram_variant_calling_tumor_only = ch_cram_variant_calling_tumor_filtered.transpose().map{ it -> [it[1], it[2], it[3]] }

            if(params.only_paired_variant_calling){
                // Normal only samples

                // 1. Join with tumor samples, in each channel there is one key per patient now. Patients without matched tumor end up with: [patient1, [meta1], [cram1,crai1], null] as there is only one matched normal possible
                ch_cram_variant_calling_normal_joined = ch_cram_variant_calling_normal_to_cross.join(ch_cram_variant_calling_tumor_grouped, remainder: true)

                // 2. Filter out entries with last entry null
                ch_cram_variant_calling_normal_filtered = ch_cram_variant_calling_normal_joined.filter{ it ->  !(it.last()) }

                // 3. Remove patient ID field & null value for further processing [meta1, [cram1,crai1]] [meta2, [cram2,crai2]] (no transposing needed since only one normal per patient ID)
                ch_cram_variant_calling_status_normal = ch_cram_variant_calling_normal_filtered.map{ it -> [it[1], it[2], it[3]] }

            }else{
                ch_cram_variant_calling_status_normal = ch_cram_variant_calling_status.normal
            }

            // Tumor - normal pairs
            // Use cross to combine normal with all tumor samples, i.e. multi tumor samples from recurrences
            ch_cram_variant_calling_pair = ch_cram_variant_calling_normal_to_cross.cross(ch_cram_variant_calling_pair_to_cross)
                .map { normal, tumor ->
                    def meta = [:]
                    meta.patient    = normal[0]
                    meta.normal_id  = normal[1].sample
                    meta.tumor_id   = tumor[1].sample
                    meta.sex        = normal[1].sex
                    meta.id         = "${meta.tumor_id}_vs_${meta.normal_id}".toString()

                    [meta, normal[2], normal[3], tumor[2], tumor[3]]
                }

            // GERMLINE VARIANT CALLING
            // TODO

            // TUMOR ONLY VARIANT CALLING
            // TODO
            
            // PAIR VARIANT CALLING
            PAIR_VARIANT_CALLING(
                params.tools,
                ch_cram_variant_calling_pair,
                dbsnp,
                dbsnp_tbi,
                dict,
                fasta,
                fasta_fai,
                germline_resource,
                germline_resource_tbi,
                intervals,
                intervals_bed_gz_tbi,
                intervals_bed_combined,
                msisensorpro_scan,
                pon,
                pon_tbi
            )

            // Gather vcf files for annotation and QC
            vcf_to_annotate = Channel.empty()
            vcf_to_annotate = vcf_to_annotate.mix(PAIR_VARIANT_CALLING.out.mutect2_vcf)

            // Gather used softwares versions
            ch_versions = ch_versions.mix(PAIR_VARIANT_CALLING.out.versions)

            // VCF QC
            VCF_QC(vcf_to_annotate, intervals_bed_combined)

            ch_versions = ch_versions.mix(VCF_QC.out.versions)
            ch_reports  = ch_reports.mix(VCF_QC.out.bcftools_stats.collect{meta, stats -> stats})
            ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_counts.collect{ meta, counts -> counts})
            ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_qual.collect{ meta, qual -> qual })
            ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_filter_summary.collect{meta, summary -> summary})

            VARIANTCALLING_CSV(vcf_to_annotate)
        }

        // ANNOTATE
        if (params.step == 'annotate') vcf_to_annotate = ch_input_sample

        if (params.tools.split(',').contains('merge') || params.tools.split(',').contains('snpeff') || params.tools.split(',').contains('vep')) {

            vep_fasta = (params.vep_include_fasta) ? fasta : []

            ANNOTATE(
                vcf_to_annotate,
                vep_fasta,
                params.tools,
                snpeff_db,
                snpeff_cache,
                vep_genome,
                vep_species,
                vep_cache_version,
                vep_cache,
                vep_extra_files)

            // Gather used softwares versions
            ch_versions = ch_versions.mix(ANNOTATE.out.versions)
            ch_reports  = ch_reports.mix(ANNOTATE.out.reports)
        }
    }



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines."
            System.exit(1)
        }
    }

    // Additional check of sample sheet:
    // 1. If params.step == "mapping", then each row should specify a lane and the same combination of patient, sample and lane shouldn't be present in different rows.
    // 2. The same sample shouldn't be listed for different patients.
    def patient_sample_lane_combinations_in_samplesheet = []
    def sample2patient = [:]

    Channel.from(csv_file).splitCsv(header: true)
        .map{ row ->
            if (!(row.patient && row.sample)){
                log.error "Missing field in csv file header. The csv file must have fields named 'patient' and 'sample'."
                System.exit(1)
            }
            if (params.step == "mapping") {
                if ( !row.lane ) {  // This also handles the case where the lane is left as an empty string
                    log.error('The sample sheet should specify a lane for patient "' + row.patient.toString() + '" and sample "' + row.sample.toString() + '".')
                    System.exit(1)
                }
                def patient_sample_lane = [row.patient.toString(), row.sample.toString(), row.lane.toString()]
                if (patient_sample_lane in patient_sample_lane_combinations_in_samplesheet) {
                    log.error('The patient-sample-lane combination "' + row.patient.toString() + '", "' + row.sample.toString() + '", and "' + row.lane.toString() + '" is present multiple times in the sample sheet.')
                    System.exit(1)
                } else {
                    patient_sample_lane_combinations_in_samplesheet.add(patient_sample_lane)
                }
            }
            if (!sample2patient.containsKey(row.sample.toString())) {
                sample2patient[row.sample.toString()] = row.patient.toString()
            } else if (sample2patient[row.sample.toString()] != row.patient.toString()) {
                log.error('The sample "' + row.sample.toString() + '" is registered for both patient "' + row.patient.toString() + '" and "' + sample2patient[row.sample.toString()] + '" in the sample sheet.')
                System.exit(1)
            }
        }

    sample_count_all = 0
    sample_count_normal = 0
    sample_count_tumor = 0

    Channel.from(csv_file).splitCsv(header: true)
        //Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            sample_count_all++
            [[row.patient.toString(), row.sample.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing

        def meta = [:]

        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no sex specified, sex is not considered
        // sex is only mandatory for somatic CNV
        if (row.sex) meta.sex = row.sex.toString()
        else meta.sex = 'NA'

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        if (meta.status == 0) sample_count_normal++
        else sample_count_tumor++

        // Two checks for ensuring that the pipeline stops with a meaningful error message if
        // 1. the sample-sheet only contains normal-samples, but some of the requested tools require tumor-samples, and
        // 2. the sample-sheet only contains tumor-samples, but some of the requested tools require normal-samples.
        if ((sample_count_normal == sample_count_all) && params.tools) { // In this case, the sample-sheet contains no tumor-samples
            def tools_tumor = ['mutect2', 'msisensorpro']
            def tools_tumor_asked = []
            tools_tumor.each{ tool ->
                if (params.tools.split(',').contains(tool)) tools_tumor_asked.add(tool)
            }
            if (!tools_tumor_asked.isEmpty()) {
                log.error('The sample-sheet only contains normal-samples, but the following tools, which were requested with "--tools", expect at least one tumor-sample : ' + tools_tumor_asked.join(", "))
                System.exit(1)
            }
        } else if ((sample_count_tumor == sample_count_all) && params.tools) {  // In this case, the sample-sheet contains no normal/germline-samples
            def tools_requiring_normal_samples = ['ascat', 'deepvariant', 'haplotypecaller', 'msisensorpro']
            def requested_tools_requiring_normal_samples = []
            tools_requiring_normal_samples.each{ tool_requiring_normal_samples ->
                if (params.tools.split(',').contains(tool_requiring_normal_samples)) requested_tools_requiring_normal_samples.add(tool_requiring_normal_samples)
            }
            if (!requested_tools_requiring_normal_samples.isEmpty()) {
                log.error('The sample-sheet only contains tumor-samples, but the following tools, which were requested by the option "tools", expect at least one normal-sample : ' + requested_tools_requiring_normal_samples.join(", "))
                System.exit(1)
            }
        }

        // mapping with fastq
        if (row.lane && row.fastq_2) {
            meta.id         = "${row.sample}-${row.lane}".toString()
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''

            def flowcell    = flowcellLaneFromFastq(fastq_1)
            //Don't use a random element for ID, it breaks resuming
            def read_group  = "\"@RG\\tID:${flowcell}.${row.sample}.${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.patient}_${row.sample}\\tLB:${row.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = 'fastq'

            meta.size       = 1 // default number of splitted fastq

            if (params.step == 'mapping') return [meta, [fastq_1, fastq_2]]
            else {
                log.error "Samplesheet contains fastq files but step is `$params.step`. Please check your samplesheet or adjust the step parameter."
                System.exit(1)
            }

        // start from BAM
        } else if (row.bam) {
            if (!row.bai) {
                log.error "BAM index (bai) should be provided."
            }
            meta.id         = "${row.sample}-${row.lane}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def bai         = file(row.bai,   checkIfExists: true)
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.sample}_${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.sample}\\tLB:${row.sample}\\tPL:${params.seq_platform}\""

            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = 'bam'

            meta.size       = 1 // default number of splitted fastq

            if (params.step != 'annotate') return [meta, bam, bai]
            else {
                log.error "Samplesheet contains bam files but step is `annotate`. The pipeline is expecting vcf files for the annotation."
                System.exit(1)
            }

        // recalibration
        } else if (row.table && row.cram) {
            meta.id   = meta.sample
            def cram  = file(row.cram,  checkIfExists: true)
            def crai  = file(row.crai,  checkIfExists: true)
            def table = file(row.table, checkIfExists: true)

            meta.data_type  = 'cram'

            if (!(params.step == 'mapping' || params.step == 'annotate')) return [meta, cram, crai, table]
            else {
                log.error "Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
                System.exit(1)
            }

        // recalibration when skipping MarkDuplicates
        } else if (row.table && row.bam) {
            meta.id   = meta.sample
            def bam   = file(row.bam,   checkIfExists: true)
            def bai   = file(row.bai,   checkIfExists: true)
            def table = file(row.table, checkIfExists: true)

            meta.data_type  = 'bam'

            if (!(params.step == 'mapping' || params.step == 'annotate')) return [meta, bam, bai, table]
            else {
                log.error "Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
                System.exit(1)
            }

        // prepare_recalibration or variant_calling
        } else if (row.cram) {
            meta.id = meta.sample
            def cram = file(row.cram, checkIfExists: true)
            def crai = file(row.crai, checkIfExists: true)

            meta.data_type  = 'cram'

            if (!(params.step == 'mapping' || params.step == 'annotate')) return [meta, cram, crai]
            else {
                log.error "Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
                System.exit(1)
            }

        // prepare_recalibration when skipping MarkDuplicates or `--step markduplicates`
        } else if (row.bam) {
            meta.id = meta.sample
            def bam = file(row.bam, checkIfExists: true)
            def bai = file(row.bai, checkIfExists: true)

            meta.data_type  = 'bam'

            if (!(params.step == 'mapping' || params.step == 'annotate')) return [meta, bam, bai]
            else {
                log.error "Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
                System.exit(1)
            }

        // annotation
        } else if (row.vcf) {
            meta.id = meta.sample
            def vcf = file(row.vcf, checkIfExists: true)

            meta.data_type     = 'vcf'
            meta.variantcaller = row.variantcaller ?: ''

            if (params.step == 'annotate') return [meta, vcf]
            else {
                log.error "Samplesheet contains vcf files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
                System.exit(1)
            }
        } else {
            log.error "Missing or unknown field in csv file header. Please check your samplesheet"
            System.exit(1)
        }
    }
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    } else {
        fcid = "NA"
    }
    return fcid
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
