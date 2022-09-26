//
// This file holds several functions specific to the workflow/stpp.nf in the stpp pipeline
//

class WorkflowSTPP {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)

        if (!params.fasta) {
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }

    public static String retrieveInput(params, log){
        switch (params.step) {
            case 'mapping':                 log.warn "Can't start with step $params.step without samplesheet"
                                            System.exit(1);
                                            break
            case 'markduplicates':          log.warn "Using file ${params.outdir}/csv/mapped.csv"
                                            params.replace("input","${params.outdir}/csv/mapped.csv");
                                            break
            case 'prepare_recalibration':   log.warn "Using file ${params.outdir}/csv/markduplicates_no_table.csv"
                                            params.replace("input", "${params.outdir}/csv/markduplicates_no_table.csv");
                                            break
            case 'recalibrate':             log.warn "Using file ${params.outdir}/csv/markduplicates.csv"
                                            params.replace("input", "${params.outdir}/csv/markduplicates.csv");
                                            break
            case 'variant_calling':         log.warn "Using file ${params.outdir}/csv/recalibrated.csv"
                                            params.replace("input", "${params.outdir}/csv/recalibrated.csv");
                                            break
            case 'annotate':                log.warn "Using file ${params.outdir}/csv/variantcalled.csv"
                                            params.replace("input","${params.outdir}/csv/variantcalled.csv");
                                            break
            default:    log.warn "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'"
                        exit 1, "Unknown step $params.step"
        }
    }
}
