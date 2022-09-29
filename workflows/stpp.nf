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
    params.mappability,
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
    if(params.pon && params.pon.contains("/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz")){
        log.warn("The default Panel-of-Normals provided by GATK is used for Mutect2.\nIt is highly recommended to generate one from normal samples that are technical similar to the tumor ones.\nFor more information: https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-")
    }
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

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[params.genome]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
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
mappability        = params.mappability        ? Channel.fromPath(params.mappability).collect()              : Channel.value([])
pon                = params.pon                ? Channel.fromPath(params.pon).collect()                      : Channel.value([]) //PON is optional for Mutect2 (but highly recommended)

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope


// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Create samplesheets to restart from different steps
include { MAPPING_CSV                                          } from '../subworkflows/local/mapping_csv'

// Build indices if needed
include { PREPARE_GENOME                                       } from '../subworkflows/local/prepare_genome'

// Run FASTQC
include { RUN_FASTQC                                           } from '../subworkflows/nf-core/run_fastqc'

// TRIM/SPLIT FASTQ Files
include { FASTP                                                } from '../modules/fastp/main'

// Map input reads to reference genome
include { GATK4_MAPPING                                        } from '../subworkflows/nf-core/gatk4/mapping/main'

// Merge and index BAM files (optional)
include { MERGE_INDEX_BAM                                      } from '../subworkflows/nf-core/merge_index_bam'












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

    // STEP 2: 






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
