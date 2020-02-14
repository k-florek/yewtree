#!/usr/bin/env nextflow

//Description: Workflow for building SNP and Core-genome Trees from raw illumina reads
//Author: Kelsey Florek
//eMail: kelsey.florek@slh.wisc.edu

//starting parameters
params.reads = "raw_reads/*R{1,2}*.fastq.gz"
params.singleEnd = false
params.outdir = "results"

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_trimming }

//Step1a: FastQC
process fastqc {
  tag "$name"
  publishDir "${params.outdir}/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  set val(name), file(reads) from read_files_fastqc

  output:
  file "*_fastqc.{zip,html}" into fastqc_results

  script:
  """
  fastqc -q  ${reads}
  """
}

//Step1b: Trim with Seqyclean
process seqyclean {
  tag "$name"
  publishDir "${params.outdir}/trimmed", mode: 'copy'

  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}_{PE1,PE2}.fastq.gz") into trimmed_reads
  file "${name}_SummaryStatistics.txt" into seqyclean_report

  script:
  if(params.singleEnd){
    """
      seqyclean -minlen 25 -qual -c /Adapters_plus_PhiX_174.fasta -U ${reads} -gz -o ${name}
    """
  }
  else {
    """
    seqyclean -minlen 25 -qual -c /Adapters_plus_PhiX_174.fasta -1 ${reads[0]} -2 ${reads[1]} -gz -o ${name}
    """
  }
}

//Step2: Assemble trimmed reads with Shovill
process shovill {
  tag "$name"
  publishDir "${params.outdir}/assembled", mode: 'copy'

  memory '8 GB'
  ram=6

  input:
  set val(name), file(reads) from trimmed_reads

  output:
  tuple name, file("${name}.contigs.fa") into assembled_genomes

  shell:
  '''
  ram=`awk '/MemTotal/ { printf "%.0f \\n", $2/1024/1024 - 1 }' /proc/meminfo`
  shovill --cpus 0 --ram $ram  --outdir . --R1 !{reads[0]} --R2 !{reads[1]} --force
  mv contigs.fa !{name}.contigs.fa
  '''
}

//Step3a: Assembly Quality Report
process quast {
  tag "$name"
  publishDir "${params.outdir}/quast",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes

  output:
  file "${name}.report.txt" into quast_report

  script:
  """
  quast.py ${assembly} -o .
  mv report.txt ${name}.report.txt
  """
}

/*
//Step3b: Annotate with prokka
process prokka {
  tag "$prefix"
  publishDir "${params.outdir}/prokka",mode:'copy'

  input:
  file(assembly) "${name}.contigs.fa" from assembled_genomes

  output:
  file "*.gff" into annotated_genomes

  script:
  """
  prokka --cpu 0 --compliant --mincontiglen 500 --outdir . ${assembly}
  """
}
*/
/*
process multiqc {
  tag "$prefix"
  publishDir "${params.outdir}/MultiQC", mode: 'copy'
  echo true

  input:
  //file multiqc_config
  file (fastqc:'fastqc/*') from fastqc_results.collect()
  file ('trimmed/*') from seqyclean_report.collect()
  file ('quast/*') from quast_report.collect()


  output:
  file "*multiqc_report.html" into multiqc_report
  file "*_data"
  val prefix into multiqc_prefix

  script:
  prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
  rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
  rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
  """
  multiqc -f $rtitle $rfilename --config $multiqc_config . 2>&1
  """
}
*/
