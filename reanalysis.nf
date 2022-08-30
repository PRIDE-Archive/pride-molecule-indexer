#!/usr/bin/env nextflow
/*
========================================================================================
                         bigbio/pride-molecules-indexer
========================================================================================
 bigbio/pride-molecules-indexer Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/bigbio/pride-molecules-indexer
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 1

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run reanalysis.nf -c nextflow.config --project_accession "PXD029360" -profile conda

    Main arguments:
      --project_accession           Project accession to convert the identifications to json files.
      --reanalysis_accession        Reanalysis accession for the specific renalysis.
      --outdir                      Output directory containing the information (json) of the project
      --reanalysis_folder           Folder with the reanalysis data files

      Advanced Options:
      --pride_production_folder     (Optional) Local folder that contains the submitted files for PRIDE

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Stage config files
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

// Validate input
if (isCollectionOrArray(params.project_accession)){
  tocheck = params.project_accession[0]
} else {
  tocheck = params.project_accession
}

params.project_accession = params.project_accession ?: { log.error "No project accession provided. Make sure you have used the '--project_accession' option."; exit 1 }()
params.reanalysis_accession = params.reanalysis_accession ?: { log.error "No project reanalysis accession provided. Make sure you have used the '--reanalysis_accession' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()
params.reanalysis_folder = params.reanalysis_folder ?: { log.error "No project reanalysis folder provided. Make sure you have used the '--reanalysis_folder' option."; exit 1 }()

mztab_files = Channel.fromPath("${params.reanalysis_folder}/*.mztab")
mzid_files  = Channel.fromPath("${params.reanalysis_folder}/*.mzid")
mzid_files2  = Channel.fromPath("${params.reanalysis_folder}/*.mzid")
mzid_files2.view()
identification_files = mztab_files.concat(mzid_files)

mzmls_files = Channel.fromPath("${params.reanalysis_folder}/*.mzML", checkIfExists: true)
sdrf_files  = Channel.fromPath("${params.reanalysis_folder}/*.sdrf.tsv", checkIfExists: true)

process generate_json_index_files{

  label 'process_high'

  publishDir "${params.outdir}", mode: 'copy', pattern: '**.json'

  input:
    file(input_files) from mzmls_files.collect()
    file(id_file) from identification_files
    file(sdrf) from sdrf_files

  output:
    file("**.json") optional true into final_index_json

  script:
  java_mem = "-Xmx" + task.memory.toGiga() + "G"
  """
  java $java_mem -jar ${baseDir}/bin/pride-molecules-indexer-1.0.0-SNAPSHOT-bin.jar generate-index-files --app.result-file="${id_file}" --app.folder-output=`pwd` --app.spectra-files="${input_files.join(",")}" --app.project-accession=${params.project_accession} --app.sample-file=${sdrf} --app.reanalysis-accession=${params.reanalysis_accession}
  """
}

// final_index_json.subscribe { println "value: $it" println " "}

//--------------------------------------------------------------- //
//---------------------- Nextflow specifics --------------------- //
//--------------------------------------------------------------- //


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['User']             = workflow.userName
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'pride-molecules-indexer-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'bigbio/pride-molecules-indexer Workflow Summary'
    section_href: 'https://github.com/bigbio/pride-molecules-indexer'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {

    label 'downloading_thread'

    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    echo $workflow.manifest.version &> v_msstats_plfq.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {

    label 'downloading_thread'

    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[bigbio/pride-molecules-indexer] Successful: $workflow.runName ${params.project_accession}"
    if (!workflow.success) {
        subject = "[bigbio/pride-molecules-indexer] FAILED: $workflow.runName ${params.project_accession}"
    }

    def msg = """\
            Pipeline execution summary: bigbio/pride-molecules-indexer
            ---------------------------
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            Error message: ${workflow.errorMessage}
            """
            .stripIndent()

    sendMail(to: params.email, subject: subject, body: msg)
}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  bigbio/pride-molecules-indexer v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}


//--------------------------------------------------------------- //
//---------------------- Utility functions  --------------------- //
//--------------------------------------------------------------- //

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Check class of an Object for "List" type
boolean isCollectionOrArray(object) {
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}