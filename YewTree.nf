#!/usr/bin/env nextflow

//Description: Workflow for building SNP and Core-genome Trees from raw illumina reads
//Author: Kelsey Florek
//eMail: kelsey.florek@slh.wisc.edu

//starting parameters
params.reads = "raw_reads/*{R1,R2,_1,_2}*.fastq.gz"
params.singleEnd = false
params.outdir = "results"

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { raw_reads }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(oldName), file(reads) from raw_reads

  output:
  tuple name, file("${name}_{R1,R2}.fastq.gz") into read_files_fastqc, read_files_trimming

  script:
  name = oldName.split("_")[0]
  if(params.singleEnd){
  """
  mv ${reads} ${name}.fastq.gz
  """
  }
  else {
  """
  mv ${reads[0]} ${name}_R1.fastq.gz
  mv ${reads[1]} ${name}_R2.fastq.gz
  """
  }
}

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

//Step1b: Trim with Trimmomatic
process trim {
  tag "$name"
  publishDir "${params.outdir}/trimmed", mode: 'copy'

  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}*{_1,_2}.fastq.gz") into trimmed_reads

  script:
  //trimming parameters
  minlength=75
  windowsize=4
  qualitytrimscore=30
  threads=4

  //TODO need to setup output for single end stuff
  if(params.singleEnd){
    """
    java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads ${threads} ${reads} -baseout ${name}.fastq.gz SLIDINGWINDOW:${windowsize}:${qualitytrimscore} MINLEN:${minlength}
    """
  }
  else {
    """
    java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${threads} ${reads} -baseout ${name}.fastq.gz SLIDINGWINDOW:${windowsize}:${qualitytrimscore} MINLEN:${minlength}
    mv ${name}*1P.fastq.gz ${name}_1.fastq.gz
    mv ${name}*2P.fastq.gz ${name}_2.fastq.gz
    """
  }
}
//Step2: Remove PhiX contamination
process cleanreads {
  tag "$name"
  publishDir "${params.outdir}/cleanedreads", mode: 'copy'

  input:
  set val(name), file(reads) from trimmed_reads

  output:
  tuple name, file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads
  file("${name}.{phix,adapters}.stats.txt") into read_cleanning_stats

  script:
  """
  repair.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_1.clean.fastq.gz out2=${name}_2.clean.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  """
}

//Step3: Assemble trimmed reads with Shovill
process shovill {
  tag "$name"
  publishDir "${params.outdir}/assembled", mode: 'copy'

  input:
  set val(name), file(reads) from cleaned_reads

  output:
  tuple name, file("${name}.contigs.fa") into assembled_genomes_quality, assembled_genomes_annotation

  shell:
  '''
  ram=`awk '/MemTotal/ { printf "%.0f \\n", $2/1024/1024 - 1 }' /proc/meminfo`
  shovill --cpus 0 --ram $ram  --outdir . --R1 !{reads[0]} --R2 !{reads[1]} --force
  mv contigs.fa !{name}.contigs.fa
  '''
}

//Step4a: Assembly Quality Report
process quast {
  publishDir "${params.outdir}/quast",mode:'copy'

  input:
  file(assemblies) from assembled_genomes_quality.collect()

  output:
  file "assembly.report.txt" into quast_report

  script:
  """
  quast.py ${assemblies} -o .
  mv report.txt assembly.report.txt
  """
}


//Step4b: Annotate with prokka
process prokka {
  tag "$name"
  publishDir "${params.outdir}/annotated",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_annotation

  output:
  file "${name}.gff" into annotated_genomes

  script:
  """
  prokka --cpu 0 --force --compliant --prefix ${name} --mincontiglen 500 --outdir . ${assembly}
  """
}

//Step5: Align with Roary
process roary {
  publishDir "${params.outdir}/core_alignment",mode:'copy'

  input:
  file(annotatedGenomes) from annotated_genomes.collect()

  output:
  file "core_gene_alignment.aln" into core_aligned_genomes
  file "roary_alignment_statistics.txt" into core_aligned_stats

  shell:
  '''
  cpus=`grep -c ^processor /proc/cpuinfo`
  roary -e -p $cpus *.gff
  mv summary_statistics.txt roary_alignment_statistics.txt
  '''
}

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
