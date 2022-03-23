//#! $HOME/miniconda3/envs/rMAP-1.0/bin/nextflow

/*
========================================================================================
               R M A P   P I P E L I N E
========================================================================================
              NEXTFLOW PIPELINE FOR rMAP
 
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    ============================================================
    main  ~  version ${params.version}
    ============================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main nf-rMAP.nf <args>
    
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --profile                     Hardware config to use. (local.config)
      --reference                   Path to reference genome to be used
      --output                      Path to output data (must be surrounded with quotes)  

    """.stripIndent()
}

// Show help emssage
params.help = false
if (params.help){
    helpMessage();
    exit 0
}

// Configurable variables
params.input = "./sequences/*.gz"
read_data = Channel.fromPath(params.input)
params.clean_reads = './results/trimmed_reads/*_{1,2}.clean.fastq.gz'
read_data = Channel.fromPath(params.input)
params.outdir = "./results"
params.reads = "./sequences/*_{1,2}.fastq.gz"
params.adapter = "/${HOME}/miniconda3/envs/rMAP-1.0/config-files/adapters.fa"
params.contigs = "./results/assembly/*.fa"
params.max_cpus = 12
params.quality = 27
params.minlength = 80
params.version = 1.0


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Header log info
log.info "============================================================"
log.info "************************************************************"
log.info " "
log.info " nextflow-rMAP pipeline  ~  version ${params.version}"
log.info " "
log.info "************************************************************"
log.info "============================================================"
def summary = [:]
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Max CPUs']       = params.max_cpus
summary['Output dir']     = params.outdir
summary['Q30 Stats']      = params.quality
summary['MinSeq Length']  = params.minlength
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile

log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "============================================================"

Channel.fromFilePairs( params.reads )
        .into { fastqc_input_ch; trimming_input_ch }

process generateList {
    publishDir "${params.outdir}", mode: "copy", overwrite: false
    
    input:
    file x from read_data
    
    output:
    file('temp.txt')
    file('list.txt')


    """
    ls ${x} >>temp.txt
    cat temp.txt | cut -d_ -f1 | sort | uniq >list.txt
    """
 }

// Run FastQC on the raw sequences
process runFastQC {
    tag { "runFastqc.${pairId}" }
    publishDir "${params.outdir}/quality_reports", mode: "copy", overwrite: false

    input:
    set pairId, file(in_fastq) from fastqc_input_ch

    output:
    file("${pairId}_fastqc/*.zip") into fastqc_output_ch

    """
    mkdir ${pairId}_fastqc
    fastqc --outdir ${pairId}_fastqc \
    ${in_fastq.get(0)} \
    ${in_fastq.get(1)}
    """
}

// Run multiqc to aggregate the fastqc output
process runMultiQC{
    publishDir "${params.outdir}/quality_reports", mode: "copy", overwrite: false

    input:
    file('*') from fastqc_output_ch.collect()

    output:
    file('multiqc_report.html')

    """
    multiqc .
    """
}

// Trim-off illumina adapters
process adapterTrimming {
  publishDir "${params.outdir}/trimmed_reads", mode: "copy", overwrite: false
  cpus 8
  memory '8 GB'

  input:
  set sampleid, file(in_fastq) from trimming_input_ch

  output:
 file"${sampleid}_1.clean.fastq.gz" into processed_reads1_ch
 file"${sampleid}_2.clean.fastq.gz" into processed_reads2_ch

  """
  trimmomatic PE -threads ${params.max_cpus} -phred33 ${in_fastq.get(0)} ${in_fastq.get(1)} ${sampleid}_1.clean.fastq.gz /dev/null ${sampleid}_2.clean.fastq.gz /dev/null ILLUMINACLIP:${params.adapter}:1:30:11 LEADING:${params.quality} TRAILING:${params.quality} MINLEN:${params.minlength}
  """
}

Channel
    .fromFilePairs( params.clean_reads )
    .set { clean_reads_ch }
// De novo assembly using shovill
process denovoAssembly {
  publishDir "${params.outdir}/assembly/", mode: "copy", overwrite: false
  cpus 8
  memory '8 GB'

  input:
  set pairId, file(in_fastq) from clean_reads_ch

  output:

file("*") into assembly_ch

  """
    mkdir -p ${pairId}
    shovill --R1 ${in_fastq.get(0)} --R2 ${in_fastq.get(1)} --cpus ${params.max_cpus} --gsize 3.4M --ram 8 --force --outdir ${pairId}
    mv ${pairId}/contigs.fa ${pairId}/${pairId}.fa
    cp ${pairId}/${pairId}.fa .
  """

}

datasets = Channel
                  .fromPath(params.contigs)

// Works but requires 'annotation.sampleid.fa' workup
process genomeAnnotation {
  publishDir "${params.outdir}/annotation/${contig_fasta}", mode: "copy", overwrite: false
  
  input:
  file contig_fasta from datasets

  output:
  file("*") into annotated_ch

  script:
  """
  prokka --outdir . --prefix ${contig_fasta} --cpus 12 ${contig_fasta} --locustag ${contig_fasta} --force
  
  """
}
