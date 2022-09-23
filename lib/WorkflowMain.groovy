//
// This file holds several functions specific to the main.nf workflow in the stpp pipeline
//

class WorkflowMain {

	//
    // Print help to screen if required
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh38 -profile singularity"
        def help_string = ''
        help_string += command
        // help_string += Schema.paramsHelp(workflow, params, command)
        // help_string += Template.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            Schema.validateParameters(workflow, params, log)
        }

        // Check input has been provided
        if (!params.input) {
            log.warn "No samplesheet specified, attempting to restart with csv files"
        }
    }
}