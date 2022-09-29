//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWA_MEM } from '../../../../modules/bwa/mem/main'

workflow GATK4_MAPPING {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ch_map_index // channel: [mandatory] mapping index
        sort         // boolean: [mandatory] true -> sort, false -> don't sort

    main:

    // ch_bam_mapped = Channel.empty()
    // ch_versions = Channel.empty()
    // ch_reports  = Channel.empty()

    // Only one of the following should be run
    BWA_MEM(ch_reads,   ch_map_index, sort) // If aligner is bwa-mem

    // ch_bam_mapped = ch_bam_mapped.mix(BWA_MEM.out.bam)
    // ch_reports = ch_reports.mix(BWA_MEM.out.log)
    // ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    emit:
        // bam      = ch_bam_mapped // channel: [ [meta], bam ]
        // reports  = ch_reports
        // versions = ch_versions   // channel: [ versions.yml ]
        bam      = BWA_MEM.out.bam // channel: [ [meta], bam ]
        reports  = BWA_MEM.out.log
        versions = BWA_MEM.out.versions.first()   // channel: [ versions.yml ]
}
