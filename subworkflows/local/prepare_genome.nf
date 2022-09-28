//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { BWA_INDEX                              } from '../../modules/bwa/index/main'
include { GATK4_CREATESEQUENCEDICTIONARY         } from '../../modules/gatk4/createsequencedictionary/main'
include { MSISENSORPRO_SCAN                      } from '../../modules/msisensorpro/scan/main'
include { SAMTOOLS_FAIDX                         } from '../../modules/samtools/faidx/main'
include { TABIX_TABIX as TABIX_DBSNP             } from '../../modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE } from '../../modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS      } from '../../modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_SNPS        } from '../../modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PON               } from '../../modules/tabix/tabix/main'


workflow PREPARE_GENOME {
    take:
        dbsnp                   // channel: [optional]  dbsnp
        fasta                   // channel: [mandatory] fasta
        fasta_fai               // channel: [optional]  fasta_fai
        germline_resource       // channel: [optional]  germline_resource
        known_indels            // channel: [optional]  known_indels
        known_snps              // channel: [optional]  known_snps
        pon                     // channel: [optional]  pon

    main:

    ch_versions = Channel.empty()

    BWA_INDEX(fasta)     // If aligner is bwa-mem

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    MSISENSORPRO_SCAN(fasta.map{ it -> [[id:it[0].baseName], it] })
    SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].getName()], it] })

    // the following are flattened and mapped in case the user supplies more than one value for the param
    // written for KNOWN_INDELS, but preemptively applied to the rest
    // [file1,file2] becomes [[meta1,file1],[meta2,file2]]
    // outputs are collected to maintain a single channel for relevant TBI files
    TABIX_DBSNP(dbsnp.flatten().map{ it -> [[id:it.baseName], it] })
    TABIX_GERMLINE_RESOURCE(germline_resource.flatten().map{ it -> [[id:it.baseName], it] })
    TABIX_KNOWN_SNPS( known_snps.flatten().map{ it -> [[id:it.baseName], it] } )
    TABIX_KNOWN_INDELS( known_indels.flatten().map{ it -> [[id:it.baseName], it] } )
    TABIX_PON(pon.flatten().map{ it -> [[id:it.baseName], it] })

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(MSISENSORPRO_SCAN.out.versions)
    ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
    ch_versions = ch_versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    ch_versions = ch_versions.mix(TABIX_KNOWN_SNPS.out.versions)
    ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions)
    ch_versions = ch_versions.mix(TABIX_PON.out.versions)

    emit:
        bwa                              = BWA_INDEX.out.index                                                 // path: bwa/*
        dbsnp_tbi                        = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()             // path: dbsnb.vcf.gz.tbi
        dict                             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                             // path: genome.fasta.dict
        fasta_fai                        = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                    // path: genome.fasta.fai
        germline_resource_tbi            = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [tbi] }.collect() // path: germline_resource.vcf.gz.tbi
        known_snps_tbi                   = TABIX_KNOWN_SNPS.out.tbi.map{ meta, tbi -> [tbi] }.collect()        // path: {known_indels*}.vcf.gz.tbi
        known_indels_tbi                 = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()      // path: {known_indels*}.vcf.gz.tbi
        msisensorpro_scan                = MSISENSORPRO_SCAN.out.list.map{ meta, list -> [list] }              // path: genome_msi.list
        pon_tbi                          = TABIX_PON.out.tbi.map{ meta, tbi -> [tbi] }.collect()               // path: pon.vcf.gz.tbi
        versions                         = ch_versions                                                         // channel: [ versions.yml ]
}