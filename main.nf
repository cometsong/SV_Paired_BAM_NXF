#!/usr/bin/env nextflow

import Helpers
import Logos

logo = new Logo()
println logo.show()

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Set Parameter Defaults ~ ~ ~ ~ ~ ~
def setParamDefaults() {
    params.help = false

    // Mandatory variable parameters
    params.bam = null
    params.bai = null

    // Configurable variable parameters
    params.annotateMergedSVifExon = true

    params.threads = 8

    // Email check:
    //if (!params.email) { exit 1, "Must supply params.email address to send pipeline report."}

    /*
    SURVIVOR merge:
        - File with VCF names and paths
        - max distance between breakpoints (0-1 percent of length, 1- number of bp) 
        - Minimum number of supporting caller
        - Take the type into account (1==yes, else no)
        - Take the strands of SVs into account (1==yes, else no)
        - Estimate distance based on the size of SV (1==yes, else no).
        - Minimum size of SVs to be taken into account.
        - Output VCF filename
    */
    params.mergeMaxBreakDist = 1000
    params.mergeMinSupportCaller = 1
    //NOTE: We use `Minimum number of supporting caller =1` rather than `=0` to avoid an issue with Delly passing 'unsupported' calls.
    params.mergeTakeTypeAccount = 1
    params.mergeTakeStrandAccount = 0
    params.mergeEstDistSizeSV = 0
    params.mergeMinSizeSV = 0

    /*
    params:
    SURVIVOR vcftobed:
        - vcf file
        - min size (-1 == no limit)
        - max size (-1 == no limit)
        - output file
    */
    params.vcftoBedMinSize = -1
    params.vcftoBedMaxSize = -1

    params.keep_intermediate = false
}
setParamDefaults()

def helpMessage() {
    log.info"""
    =========================================
      ${workflow.manifest.name} v${workflow.manifest.version}
    =========================================
    ${workflow.manifest.description}

    NOTE: This will create a subdirectory *within* your current directory ($PWD)
          named using the timestamp and the sample name from the bam file you submit.

    Usage:
      The typical command for running the pipeline is as follows:
        nextflow run ${workflow.manifest.name} -profile sumner --bam "/path/to/file.bam" --bai "/path/to/file.bai"
      OR if using a config file with these parameters:
        nextflow -c your_params.config run ${workflow.manifest.name} -profile sumner
    Mandatory:
      --bam       Path to duplicate filtered BAM file.
      --bai       Path to associated BAM index file.
    Optional:
      --annotateMergedSVifExon
                  Annotates merged SV call if is in an exon. 
      --threads   Number of cpus to use for multi-threaded processes.
      --email     The email address to send the pipeline report.
      -name       Name for the pipeline run. If not specified Nextflow will
                      automatically generate a random mnemonic.
      -profile    Environment config to use.
                      [ choices: standard (local), sumner, helix, winter ]
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Summary Info ~ ~ ~ ~ ~ ~
// Header info
def summary = [:]
summary['Pipeline']         = workflow.manifest.name
summary['Description']      = workflow.manifest.description
if(workflow.revision) {
    summary['Pipeline Release'] = workflow.revision
}
summary['Run Name']         = workflow.runName
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
summary['Config Files']     = workflow.configFiles
summary['Command Line']     = workflow.commandLine
summary['Nextflow Info']    = "v${nextflow.version}, build: ${nextflow.build}, on ${nextflow.timestamp}"
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
if(workflow.containerEngine) {
    summary['Container Engine'] = "$workflow.containerEngine"
}
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus"

// Pipeline Params:
summary['Parameters......'] = ''
if(params.email) {
    summary['. E-mail']     = params.email }
summary['. BAM file']       = params.bam
summary['. BAM index']      = params.bai
summary['. Annotate if Exon'] = params.annotateMergedSVifExon

//~~~~~ Parse Name of Sample ~~~~~\\
def bam_file_name = file("${params.bam}").name
def sample_name_RE = ~/(.*?)(_|\.)*.bam/
def sample_name_bits = bam_file_name =~ sample_name_RE
if ( sample_name_bits.size() >= 1) {
    sample_name = sample_name_bits[0][1]
} else {
    sample_name = "no sample_name matched... oops"
}
summary['. Sample Name'] = sample_name

def timestamp = new Date().format("yyyyMMdd'-'hhmmss")
def outdir = "${timestamp}_${sample_name}"
summary['. Output dir'] = outdir

summary['Run Start Time']   = workflow.start
println Summary.show(summary)

//~~~~~ Get Abspath of 'outdir' ~~~~~\\
abs_outdir = file(outdir).toAbsolutePath().toString()

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Extra Reference Files ~~~~~
def SV_ref_dir = "$baseDir/refs"

//FIXME: Copy all ref files to single dir in pipeline location (portability) e.g. the bam index folder:
//def index_dir_bwa = SV_ref_dir + "/GRCm38_BWAIndex/GRCm38.p5_full.fa"
def index_dir_bwa = "/projects/lloydm/SV_sim/keshar/mmuRef/GRCm38.p5/BWAIndex/GRCm38.p5_full.fa"

def black_regions_tsv = SV_ref_dir + "/mouse.mm10.excl.tsv"
def black_regions_bed = SV_ref_dir + "/mm10.gaps.centro_telo.scafold.exclude.bed"
def mouse_ref_fasta   = SV_ref_dir + "/GRCm38_p5_primaryAssembly_rename.fa"
def mm10_exon_list    = SV_ref_dir + "/mm10_exon_list.bed"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Opening BAM/I Channel ~~~~~
Channel
    .of( sample_name, params.bam, params.bai )
    .toList()
    .into{ in_brkdncr; in_delly; in_lumpy; in_manta }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BreakDancer SV ~~~~~
process breakdancer_calling_sv {
    tag "$sample_name"
    label 'breakdancer'
    publishDir "${outdir}/BreakDancerSVOut", pattern: "*_BreakDancer-SV", mode: 'move'
    publishDir "${outdir}/BreakDancerSVOut", pattern: "*_BreakDancerSortVCF.vcf", mode: 'move'

    input:
    tuple sample_name, bam_input, bam_index from in_brkdncr
    val abs_outdir from abs_outdir

    output:
    path(breakdancer_sv_out)
    path(breakdancer_sort_vcf)
    path("vcf_path") into vcf_breakdancer

    script:
    breakdancer_config   = sample_name + "_config"
    breakdancer_sv_out   = sample_name + "_BreakDancer-SV"
    breakdancer_2_vcf    = sample_name + "_BreakDancer2VCF.vcf"
    breakdancer_sort_vcf = sample_name + "_BreakDancerSortVCF.vcf"

    log.info "Calling BreakDancer SV"
    """
    bam2cfg.pl ${bam_input} > ${breakdancer_config}
    breakdancer-max -r 5 -s 50 -h ${breakdancer_config} > ${breakdancer_sv_out}
    breakdancer2vcfHeader.py -i ${breakdancer_sv_out} -o ${breakdancer_2_vcf}
    vcfSort.sh ${breakdancer_2_vcf} ${breakdancer_sort_vcf}
    echo ${abs_outdir}/BreakDancerSVOut/${breakdancer_sort_vcf} > vcf_path # for later merging
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Delly SV ~~~~~
process delly_calling_sv {
    tag "$sample_name"
    label 'delly'

    input:
    tuple sample_name, bam_input, bam_index from in_delly
    path(index_dir_bwa) from index_dir_bwa
    path(black_regions_tsv) from black_regions_tsv

    output:
    tuple sample_name, path(delly_bcf) into delly_bcf_out

    script:
    delly_bcf = sample_name + "_Dellybcf.bcf"
    log.info "Calling Delly SV"
    """
    delly call \
          -q 40 \
          -x ${black_regions_tsv} \
          -s 500 \
          -o ${delly_bcf} \
          -g ${index_dir_bwa} ${bam_input}
    """
}
process delly_bcf2vcf_sort {
    tag "$sample_name"
    label 'bcftools'
    publishDir "${outdir}/DellySVOut", pattern: "*_dellySort.vcf", mode: 'move'

    input:
    tuple sample_name, path(delly_bcf) from delly_bcf_out
    val abs_outdir from abs_outdir

    output:
    path(delly_sort_vcf)
    path("vcf_path") into vcf_delly

    script:
    delly_vcf      = sample_name + "_DellyVCF.vcf"
    delly_sort_vcf = sample_name + "_dellySort.vcf"

    log.info "Delly bcf2vcf to Sorted VCF"
    """
    bcftools view ${delly_bcf} > ${delly_vcf}
    vcfSort.sh ${delly_vcf} ${delly_sort_vcf}
    echo ${abs_outdir}/DellySVOut/${delly_sort_vcf} > vcf_path # for later merging
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Lumpy SV ~~~~~
process lumpy_mapping {
    tag "$sample_name"
    label 'lumpy_sams'
    publishDir "${outdir}/mapped_lumpy", pattern: "*_alignBWA_lumpy.bam", mode: 'move'
    cpus params.threads

    input:
    tuple sample_name, bam_input, bam_index from in_lumpy

    output:
    tuple sample_name, path(bam_bwa_lumpy) into ( bam_bwa_lumpy_ch, bam_bwa_lumpy_splits_ch )
    tuple sample_name, path(dis_unsorted_bam) into dis_unsorted_bam_ch

    script:
    log.info "Mapping Lumpy (Map clipped reads, read group info, extract discordant alignments)"

    // Lumpy File Names:
    bam_name_sort      = sample_name + "_alignBWA_ReadNameSort"
    bam_name_sort_full = sample_name + "_alignBWA_ReadNameSort.bam"
    bam_bwa_lumpy      = sample_name + "_alignBWA_lumpy.bam"
    bam_bwa_lumpy_sort = sample_name + "_alignBWA_lumpySort.bam"
    dis_unsorted_bam   = sample_name + "_discordants.unsorted.bam"
    dis_sorted_bam     = sample_name + "_discordants.sorted.bam"
    split_unsorted_bam = sample_name + "_splitters.unsorted.bam"
    split_sorted_bam   = sample_name + "_splitters.sorted.bam"
    """
    # Clipped_rc reads mapping to Genome
    samtools sort -n ${bam_input} -o ${bam_name_sort_full} -@ ${params.threads}
    # manual read group info
    samtools view -h ${bam_name_sort_full} \
    | samblaster --acceptDupMarks --excludeDups --addMateTags \
                 --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -@ ${params.threads} -S -b - > ${bam_bwa_lumpy}
    # Extract the discordant pairedend alignments
    samtools view -@ ${params.threads} -b -F 1294 ${bam_bwa_lumpy} > ${dis_unsorted_bam}
    """
}
process lumpy_bwa_sort {
    tag "$sample_name"
    label 'picard'

    input:
    tuple sample_name, path(bam_bwa_lumpy) from bam_bwa_lumpy_ch

    output:
    tuple sample_name, path("${bam_bwa_lumpy_sort}.ba[mi]") into bam_bwa_lumpy_sort_ch

    script:
    log.info "Mapping Lumpy (bam sort)"
    bam_bwa_lumpy_sort = sample_name + "_alignBWA_lumpySort" //+ .bam
    """
    picard SortSam I=${bam_bwa_lumpy} O=${bam_bwa_lumpy_sort}.bam \
        SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
    """
}
process lumpy_discordant_sort {
    tag "$sample_name"
    label 'picard'

    input:
    tuple sample_name, path(dis_unsorted_bam) from dis_unsorted_bam_ch

    output:
    tuple sample_name, path("${dis_sorted_bam}.ba[mi]") into dis_sorted_bam_ch

    script:
    log.info "Mapping Lumpy (discordant sort)"
    dis_sorted_bam = sample_name + "_discordants.sorted" //+ .bam
    """
    picard SortSam I=${dis_unsorted_bam} O=${dis_sorted_bam}.bam \
        SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
    """
}
process lumpy_extract_splits {
    tag "$sample_name"
    label 'lumpy_sams'

    input:
    tuple sample_name, path(bam_bwa_lumpy) from bam_bwa_lumpy_splits_ch

    output:
    tuple sample_name, path("${split_unsorted_bam}.ba[mi]") into split_unsorted_bam_ch

    script:
    log.info "Mapping Lumpy (extract the splitread alignments)"
    extractSplitReads_BwaMem = "/lumpy-sv/scripts/extractSplitReads_BwaMem"
    split_unsorted_bam = sample_name + "_splitters.unsorted" //+ .bam
    """
    samtools view -h ${bam_bwa_lumpy} \
    | ${extractSplitReads_BwaMem} -i stdin \
    | samtools view -Sb - > ${split_unsorted_bam}.bam
    """
}
process lumpy_split_bam_sort {
    tag "$sample_name"
    label 'picard'
    publishDir "${outdir}/mapped_lumpy", pattern: "*_splitters.sorted.ba*", mode: 'move'

    input:
    tuple sample_name, path(split_unsorted_bam) from split_unsorted_bam_ch

    output:
    tuple sample_name, path("${split_sorted_bam}.ba[mi]") into split_sorted_bam_ch

    script:
    log.info "Mapping Lumpy (sort split reads)"
    split_sorted_bam = sample_name + "_splitters.sorted" //+ .bam
    """
    picard SortSam I=${split_unsorted_bam[0]} O=${split_sorted_bam}.bam \
        SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
    """
}
process lumpy_call_sv {
    tag "$sample_name"
    label 'lumpy_sams'
    publishDir "${outdir}/lumpySVOut", pattern: "*_lumpySort.vcf", mode: 'move'
    publishDir "${outdir}/mapped_lumpy", pattern: "*_discordants.sorted.ba*", mode: 'move'

    //errorStrategy { task.exitStatus=141 ? 'ignore' : 'terminate' } // validExitStatus 141 for pairend_distro

    input:
    tuple sample_name, path(bam_bwa_lumpy_sort) from bam_bwa_lumpy_sort_ch
    tuple sample_name, path(split_sorted_bam) from split_sorted_bam_ch
    tuple sample_name, path(dis_sorted_bam) from dis_sorted_bam_ch
    path(black_regions_bed) from black_regions_bed
    val abs_outdir from abs_outdir

    output:
    path(dis_sorted_bam, includeInputs=true)
    path(lumpy_sort_vcf)
    path("vcf_path") into vcf_lumpy

    shell:
    log.info "Call SV by Lumpy, sort vcf"
    pairend_distro = "/lumpy-sv/scripts/pairend_distro.py"
    histo          = sample_name + "_alignBWA_lumpySort.lib1.histo"
    lumpy_vcf      = sample_name + "_lumpyOut.vcf"
    lumpy_sort_vcf = sample_name + "_lumpySort.vcf"
    '''
      RG_ID=$(samtools view -H !{bam_bwa_lumpy_sort[1]} | grep '^@RG' | sed "s/.*ID:\\([^\\t]*\\).*/\\1/g")
    #orig: metrics=$(samtools view -r "${RG_ID}" !{bam_bwa_lumpy_sort[1]} | tail -n+100000 | !{pairend_distro} -r 150 -X 4 -N 10000 -o !{histo}) 2>&1
    samtools view -r "${RG_ID}" !{bam_bwa_lumpy_sort[1]} | tail -n+100000 > pre_metrics 2>/dev/null
    metrics=$(cat pre_metrics | !{pairend_distro} -r 150 -X 4 -N 10000 -o !{histo}) 2>/dev/null \
        && [ $? = 141 ] && echo 'metrics to pairend_distro had exitcode: '$?;
       mean=$(echo "${metrics}" | cut -d " " -f 1)
       mean=$(echo "${mean}"    | cut -d ":" -f 2)
    std_dev=$(echo "${metrics}" | cut -d " " -f 2)
    std_dev=$(echo "${std_dev}" | cut -d ":" -f 2)
    rm pre_metrics;

    lumpy \
        -mw 4 \
        -x !{black_regions_bed} \
        -pe id:"${RG_ID}",bam_file:!{dis_sorted_bam[1]},histo_file:!{histo},mean:"${mean}",stdev:"${std_dev}",read_length:150,min_non_overlap:150,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
        -sr id:"${RG_ID}",bam_file:!{split_sorted_bam[1]},back_distance:10,weight:1,min_mapping_threshold:20 \
        > !{lumpy_vcf}

    vcfSort.sh !{lumpy_vcf} !{lumpy_sort_vcf}
    echo !{abs_outdir}/lumpySVOut/!{lumpy_sort_vcf} > vcf_path # for later merging
    '''
    //FIXME: check/compare results of vcfSort.sh to lumpySort.sh/BDSort/etc
    //orig: >  grep '^#' ${lumpy_vcf} > chr_vcf && grep '^chr' ${lumpy_vcf} >> chr_vcf
    //orig: >  grep '^#' chr_vcf > ${sort_vcf} && grep -v '^#' chr_vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> ${sort_vcf} 
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Manta SV ~~~~~
process manta_calling_sv {
    tag "$sample_name"
    label 'manta'
    publishDir "${outdir}/mantaSVout", pattern: "*candidateSV.vcf*", mode: 'move', \
        saveAs: { vcf -> "${sample_name}_manta_${vcf}" }
    publishDir "${outdir}/temps", pattern: "mantaSVOut", enabled: params.keep_intermediate
    cpus params.threads

    input:
    tuple sample_name, bam_input, bam_index from in_manta
    path(mouse_ref_fasta) from mouse_ref_fasta
    val abs_outdir from abs_outdir

    output:
    path("*candidateSV.vcf*")
    path("vcf_path") into vcf_manta

    script:
    log.info "Calling Manta SV"
    """
    /manta/bin/configManta.py \
        --runDir mantaSVOut \
        --bam ${bam_input} \
        --referenceFasta ${mouse_ref_fasta}
    ./mantaSVOut/runWorkflow.py -m local -j ${params.threads}
    mv mantaSVOut/results/variants/candidateSV.vcf.gz ./
    gunzip candidateSV.vcf.gz
    //TODO: use this instead of saveAs: gunzip -c candidateSV.vcf.gz > ${sample_name}_manta_candidateSV.vcf
    echo ${abs_outdir}/mantaSVout/${sample_name}_manta_candidateSV.vcf > vcf_path # for later merging
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Collect all VCF file Paths  ~~~~~
// NOTE: merged extension contains 'BDLM' => Breakdancer + Delly + Lumpy + Manta
vcf_breakdancer.concat(
vcf_delly,
vcf_lumpy,
vcf_manta
)
.collectFile(name: "sample.vcfs.txt", sort: true)
.set { sample_vcfs_paths }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Merge all SV calls VCF files  ~~~~~
process survivor_merge_sv_vcf {
    tag "$sample_name"
    label 'survivor'
    label 'survivor'
    publishDir "${outdir}", pattern: "*_mergedCall.BDLM.vcf", mode: "copy"
    publishDir "${outdir}/temps", pattern: "*.vcfs.txt", enabled: params.keep_intermediate

    input:
    path(vcf_paths) from sample_vcfs_paths

    output:
    tuple sample_name, path(outMerged) into ( vcf_merged, vcf_mrg_annot )

    script:
    log.info "Merging SV Call VCFs with SURVIVOR"

    maxBreak   = params.mergeMaxBreakDist
    minSupp    = params.mergeMinSupportCaller
    takeType   = params.mergeTakeTypeAccount
    takeStrand = params.mergeTakeStrandAccount
    estDist    = params.mergeEstDistSizeSV
    minSizeSV  = params.mergeMinSizeSV
    outMerged  = sample_name + "_mergedCall.BDLM.vcf"
    """
    SURVIVOR merge  \
        $vcf_paths  \
        $maxBreak   \
        $minSupp    \
        $takeType   \
        $takeStrand \
        $estDist    \
        $minSizeSV  \
        $outMerged;
    """
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Annotate if in Exon  ~~~~~
process survivor_vcf_to_bed {
    tag "$sample_name"
    label 'survivor'
    publishDir "${outdir}/temps", pattern: "*.bedpe", enabled: params.keep_intermediate

    when:
    params.annotateMergedSVifExon == true

    input:
    tuple sample_name, path(in_vcf) from vcf_merged

    output:
    tuple sample_name, path(out_bedpe) into vcf_bedpe

    script:
    log.info "Preparing for annotation by converting VCF to BED with SURVIVOR"
    minSizeBed = params.vcftoBedMinSize
    maxSizeBed = params.vcftoBedMaxSize
    out_bedpe  = in_vcf + ".bedpe"
    """
    SURVIVOR vcftobed \
        $in_vcf \
        $minSizeBed \
        $maxSizeBed \
        $out_bedpe;
    """
}
process annot_bedpe_to_multi {
    tag "$sample_name"
    label 'python3'
    publishDir "${outdir}/temps", pattern: "*.bedpe", enabled: params.keep_intermediate

    when:
    params.annotateMergedSVifExon == true

    input:
    tuple sample_name, path(in_bedpe) from vcf_bedpe

    output:
    tuple sample_name, path(out_bedpe) into annot_bedpe

    script:
    log.info "Convert TRA calls from single line to multi-line to facilitate bedtools intersect"
    out_bedpe = in_bedpe.getBaseName() + ".annot.bedpe"
    """
    python3 ${workflow.projectDir}/bin/bedpe_annot_para.py \
        ${in_bedpe} ${out_bedpe}
    """
}
process bedtools_intersect {
    tag "$sample_name"
    label 'bedtools'
    publishDir "${outdir}/temps", pattern: "*.bedpe", enabled: params.keep_intermediate

    when:
    params.annotateMergedSVifExon == true

    input:
    tuple sample_name, path(in_bedpe) from annot_bedpe
    path(mm10_exon_list) from mm10_exon_list

    output:
    tuple sample_name, path(out_bedpe) into intersect_bedpe

    script:
    log.info "Intersect calls against mm10_exon_list"
    out_bedpe = in_bedpe.getBaseName() + ".intersect.bedpe"
    """
    bedtools intersect \
        -a ${in_bedpe} \
        -b ${mm10_exon_list} \
        -loj > ${out_bedpe}
    """
}
process annot_merge_vcf {
    tag "$sample_name"
    label 'python3'
    publishDir "${outdir}", pattern: "*.ExonAnnot.vcf", mode: 'move'
    publishDir "${outdir}/logs", pattern: "*.error.log", mode: 'move'

    when:
    params.annotateMergedSVifExon == true

    input:
    tuple sample_name, path(in_bedpe) from intersect_bedpe
    tuple sample_name, path(merged_vcf) from vcf_mrg_annot

    output:
    tuple sample_name, path(exon_vcf) into exon_annot_vcf
    path("*.error.log") optional true into errlog_ch

    script:
    log.info "Apply 'InExon' to original VCF file. Based on matching annotations from VCF and annotated bedpe."
    exon_vcf = sample_name + '.ExonAnnot.vcf'
    """
    python3 ${workflow.projectDir}/bin/annot_vcf_with_exon.py \
        ${in_bedpe} \
        ${merged_vcf} \
        > ${exon_vcf} \
        2>${exon_vcf}.error.log
    """
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Closing Info ~ ~ ~ ~ ~ ~
workflow.onComplete {
    wfEnd = [:]
    wfEnd['Completed at'] = workflow.complete
    wfEnd['Duration']     = workflow.duration
    wfEnd['Exit status']  = workflow.exitStatus
    wfEnd['Success']      = workflow.success
    if(!workflow.success){
        wfEnd['!!Execution failed'] = ''
        wfEnd['.    Error']   = workflow.errorMessage
        wfEnd['.    Report']  = workflow.errorReport
    }
    Summary.show(wfEnd)
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// vim: set ft=groovy.nextflow ts=4 sw=0 tw=100 et fdm=syntax:
