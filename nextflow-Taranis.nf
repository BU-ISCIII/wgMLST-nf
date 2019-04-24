#!/usr/bin/env nextflow

/*
========================================================================================
                  B A C T E R I A L   W G S   P R A C T I C E
========================================================================================
 #### Homepage / Documentation
 https://github.com/BU-ISCIII/bacterial_wgs_training
 @#### Authors
 Sara Monzon <smonzon@isciii.es>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1. : Preprocessing
 	- 1.1: FastQC for raw sequencing reads quality control
 	- 1.2: Trimmomatic
 - 2. : Assembly
 	- 2.1 : Assembly with spades
 	- 2.2 : Assembly stats
 - 3. : Annotation
 	- 3.1 : Prokka automatic annotation
 - 4. : MultiQC
 - 5. : Output Description HTML
 ----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     BU-ISCIII/bacterial_wgs_training : WGS analysis practice v${version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run BU-ISCIII/bacterial_wgs_training --reads '*_R{1,2}.fastq.gz' --fasta listeria.fasta --step preprocessing
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).
      --fasta                       Path to Fasta reference
    References
      --gtf                         Path to GTF reference file. (Mandatory if step = assembly)
      --scheme                      Path to core gene scheme
    Options:
      --singleEnd                   Specifies that the input is single end reads
    Trimming options
      --notrim                      Specifying --notrim will skip the adapter trimming step.
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.
      --trimmomatic_adapters_file   Adapters index for adapter removal
      --trimmomatic_adapters_parameters Trimming parameters for adapters. <seed mismatches>:<palindrome clip threshold>:<simple clip threshold>. Default 2:30:10
      --trimmomatic_window_length   Window size. Default 4
      --trimmomatic_window_value    Window average quality requiered. Default 20
      --trimmomatic_mininum_length  Minimum length of reads
    Assembly options
    Other options:
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '0.1'

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

/*
 * Default and custom value for configurable variables
 */

params.fasta = false
if( params.fasta ){
    fasta_file = file(params.fasta)
    if( !fasta_file.exists() ) exit 1, "Fasta file not found: ${params.fasta}."
}


// gtf file
params.gtf = false

if( params.gtf ){
    gtf_file = file(params.gtf)
    if( !gtf_file.exists() ) exit 1, "GTF file not found: ${params.gtf}."
}

// Trimming default
params.notrim = false


// MultiQC config file
params.multiqc_config = "${baseDir}/conf/multiqc_config.yaml"

if (params.multiqc_config){
	multiqc_config = file(params.multiqc_config)
}

// Output md template location
output_docs = file("$baseDir/docs/output.md")

// Output files options
params.saveTrimmed = false

// Default trimming options
params.trimmomatic_adapters_file = "\$TRIMMOMATIC_PATH/adapters/NexteraPE-PE.fa"
params.trimmomatic_adapters_parameters = "2:30:10"
params.trimmomatic_window_length = "4"
params.trimmomatic_window_value = "20"
params.trimmomatic_mininum_length = "50"

// SingleEnd option
params.singleEnd = false

// Validate  mandatory inputs
params.reads = false
if (! params.reads ) exit 1, "Missing reads: $params.reads. Specify path with --reads"

if ( ! params.gtf ){
    exit 1, "GTF file not provided for assembly step, please declare it with --gtf /path/to/gtf_file"
}

/*
 * Create channel for input files
 */

// Create channel for input reads.
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { raw_reads_fastqc; raw_reads_trimming }


// Header log info
log.info "========================================="
log.info " BU-ISCIII/bacterial_wgs_training : WGS analysis practice v${version}"
log.info "========================================="
def summary = [:]
summary['Reads']               = params.reads
summary['Data Type']           = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Fasta Ref']           = params.fasta
summary['GTF File']            = params.gtf
summary['Keep Duplicates']     = params.keepduplicates
summary['Step']                = params.step
summary['Container']           = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Working dir']         = workflow.workDir
summary['Output dir']          = params.outdir
summary['Script dir']          = workflow.projectDir
summary['Save Trimmed']        = params.saveTrimmed
if( params.notrim ){
    summary['Trimming Step'] = 'Skipped'
} else {
    summary['Trimmomatic adapters file'] = params.trimmomatic_adapters_file
    summary['Trimmomatic adapters parameters'] = params.trimmomatic_adapters_parameters
    summary["Trimmomatic window length"] = params.trimmomatic_window_length
    summary["Trimmomatic window value"] = params.trimmomatic_window_value
    summary["Trimmomatic minimum length"] = params.trimmomatic_mininum_length
}
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "===================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

/*
 * STEP 1.1 - FastQC
 */
process fastqc {
	tag "$prefix"
	publishDir "${params.outdir}/fastqc", mode: 'copy',
		saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

	input:
	set val(name), file(reads) from raw_reads_fastqc

	output:
	file '*_fastqc.{zip,html}' into fastqc_results
	file '.command.out' into fastqc_stdout

	script:

	prefix = name - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(_00*)?(\.fq)?(\.fastq)?(\.gz)?$/
	"""
	fastqc -t 1 $reads
	"""
}

process trimming {
	tag "$prefix"
	publishDir "${params.outdir}/trimming", mode: 'copy',
		saveAs: {filename ->
			if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
			else if (filename.indexOf(".log") > 0) "logs/$filename"
else if (filename.indexOf(".fastq.gz") > 0) "trimmed/$filename"
			else params.saveTrimmed ? filename : null
	}

	input:
	set val(name), file(reads) from raw_reads_trimming

	output:
	file '*_paired_*.fastq.gz' into trimmed_paired_reads,trimmed_paired_reads_bwa,trimmed_paired_reads_unicycler,trimmed_paired_reads_wgsoutbreaker,trimmed_paired_reads_plasmidid,trimmed_paired_reads_mlst,trimmed_paired_reads_res,trimmed_paired_reads_sero,trimmed_paired_reads_vir
	file '*_unpaired_*.fastq.gz' into trimmed_unpaired_reads
	file '*_fastqc.{zip,html}' into trimmomatic_fastqc_reports
	file '*.log' into trimmomatic_results

	script:
	prefix = name - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(_00*)?(\.fq)?(\.fastq)?(\.gz)?$/
	"""
	trimmomatic PE -phred33 $reads -threads 1 $prefix"_paired_R1.fastq" $prefix"_unpaired_R1.fastq" $prefix"_paired_R2.fastq" $prefix"_unpaired_R2.fastq" ILLUMINACLIP:${params.trimmomatic_adapters_file}:${params.trimmomatic_adapters_parameters} SLIDINGWINDOW:${params.trimmomatic_window_length}:${params.trimmomatic_window_value} MINLEN:${params.trimmomatic_mininum_length} 2> ${name}.log
	gzip *.fastq
	fastqc -q *_paired_*.fastq.gz
	"""
}


/*
* STEP 5 Assembly
*/

process unicycler {
	tag "$prefix"
	publishDir path: { "${params.outdir}/unicycler" }, mode: 'copy'

	input:
	set file(readsR1),file(readsR2) from trimmed_paired_reads_unicycler

	output:
	file "${prefix}_assembly.fasta" into scaffold_quast,scaffold_prokka,scaffold_plasmidid,scaffold_taranis

	script:
	prefix = readsR1.toString() - ~/(.R1)?(_1)?(_R1)?(_trimmed)?(_paired)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
	"""
	unicycler -1 $readsR1 -2 $readsR2 --pilon_path \$PILON_PATH -o .
	mv assembly.fasta $prefix"_assembly.fasta"
	"""
}

process quast {
	tag "$prefix"
	publishDir path: {"${params.outdir}/quast"}, mode: 'copy',
						saveAs: { filename -> if(filename == "quast_results") "${prefix}_quast_results"}

	input:
	file scaffolds from scaffold_quast.collect()
	file fasta from fasta_file
	file gtf from gtf_file

	output:
	file "quast_results" into quast_results
	file "quast_results/latest/report.tsv" into quast_multiqc

	script:
	prefix = scaffolds[0].toString() - ~/(_scaffolds\.fasta)?$/
	"""
	quast.py -R $fasta -G $gtf $scaffolds
	"""
}

process prokka {
	tag "$prefix"
	publishDir path: {"${params.outdir}/prokka"}, mode: 'copy',
						saveAs: { filename -> if(filename == "prokka_results") "${prefix}_prokka_results"}

	input:
	file scaffold from scaffold_prokka

	output:
	file "prokka_results" into prokka_results
	file "prokka_results/${prefix}_prokka.txt" into prokka_multiqc

	script:
	prefix = scaffold.toString() - ~/(_scaffolds\.fasta)?$/
	"""
	prokka --force --outdir prokka_results --prefix prokka --genus Listeria --species monocytogenes --strain $prefix --locustag BU-ISCIII --compliant --kingdom Bacteria $scaffold
	mv prokka_results/prokka.txt prokka_results/${prefix}_prokka.txt
	"""
}


/*
 * STEP 11 MultiQC
 */


process multiqc_assembly {
	tag "$prefix"
	publishDir "${params.outdir}/MultiQC", mode: 'copy'

	input:
	file multiqc_config
	file (fastqc:'fastqc/*') from fastqc_results.collect()
	file ('trimommatic/*') from trimmomatic_results.collect()
	file ('trimommatic/*') from trimmomatic_fastqc_reports.collect()
	file ('prokka/*') from prokka_multiqc.collect()
	file ('quast/*') from quast_multiqc.collect()

	output:
	file '*multiqc_report.html' into multiqc_report
	file '*_data' into multiqc_data
	file '.command.err' into multiqc_stderr
	val prefix into multiqc_prefix

	script:
	prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'

	"""
	multiqc -d . --config $multiqc_config
	"""

}


workflow.onComplete {
	log.info "BU-ISCIII - Pipeline complete"
}



