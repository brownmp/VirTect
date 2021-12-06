version 1.0


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create the task cutadapt
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cutadapt finds and removes adapter sequences, primers,
#       poly-A tails and other types of unwanted sequence from
#       your high-throughput sequencing reads.

task RunCutadapt {
    input {
        File fastq1
        File? fastq2

        Int cpus
        Int preemptible
        String docker
        String sample_id

    }

    # Set up Prefix 
    String prefix = sample_id

    command <<<
        set -e

        # special case for tar of fastq files
        if [[ "~{fastq1}" == *.tar.gz ]] ; then
            mkdir fastq
            tar -I pigz -xvf ~{fastq1} -C fastq
            fastqs=$(find fastq -type f)
            fastq1=$fastqs[0]
            fastq2=$fastqs[1]
        fi


        python3 /usr/local/src/VirTect/VerTect_cutadapt.py \
                -1 ~{fastq1} \
                -2 ~{fastq2} \
                -F "AGATCGGAAGAG" \
                -R "AGATCGGAAGAG" \
                --cores ~{cpus}
    >>>

    output {
        File fastq1_out = "sample1_trimmed_1.fq"
        File? fastq2_out = "sample2_trimmed_2.fq"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fastq1, "GB")*3 ) + " HDD"
        docker: docker
        cpu: 1
        memory: "10GB"
    }
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create the task RunVirTect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
task RunTopHat {
    input {
        File fastq1
        File? fastq2
        
        File Virus_Reference
        File Human_Reference
        File GTF_Reference

        Int cpus
        Int preemptible
        String docker
        String sample_id
    }

    #~~~~~~~~~~~~~~~~~
    # Set the Prefix 
    #~~~~~~~~~~~~~~~~~
    String prefix = sample_id

    command <<<

        set -e

        # special case for tar of fastq files
        if [[ "~{fastq1}" == *.tar.gz ]] ; then
            mkdir fastq
            tar -I pigz -xvf ~{fastq1} -C fastq
            fastqs=$(find fastq -type f)
            fastq1=$fastqs[0]
            fastq2=$fastqs[1]
        fi
        
        
        # Untar the references  
        tar -xvf ~{Human_Reference}

        #~~~~~~~~~~~~~~~
        # TopHat
        #~~~~~~~~~~~~~~~
        tophat2 \
            -o . \
            -p ~{cpus} \
            -G ~{GTF_Reference} \
            GRCh38.genome \
            ~{fastq1} ~{fastq2}


    >>>

    output {
        File unmapped_bam = "OUTPUT/unmapped.bam"
        File accepted_hits_bam = "OUTPUT/accepted_hits.bam"
        File junctions_bed = "OUTPUT/junctions.bed"
        File deletions_bed = "OUTPUT/deletions.bed"
        File insertions_bed = "OUTPUT/insertions.bed"
        File align_summary_txt = "OUTPUT/align_summary.txt"
        File prep_reads_info = "OUTPUT/prep_reads.info"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fastq1, "GB")*3 ) + " HDD"
        docker: docker
        cpu: cpus
        memory: "10GB"
    }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create the task bam2fastq
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
task bam2fastq {
    input {
        
        File Virus_Reference
        File Human_Reference
        File GTF_Reference

        File unmapped_bam

        Int cpus
        Int preemptible
        String docker
        String sample_id
    }

    #~~~~~~~~~~~~~~~~~
    # Set the Prefix 
    #~~~~~~~~~~~~~~~~~
    String prefix = sample_id

    command <<<

        set -e

        #~~~~~~~~~~~~~~~
        # SORT
        #~~~~~~~~~~~~~~~
        samtools sort \
            -n ~{unmapped_bam} \
            -o unmapped_sorted.bam

        #~~~~~~~~~~~~~~~
        # Bedtools 
        #~~~~~~~~~~~~~~~

        bedtools bamtofastq -i . \
                unmapped_sorted.bam \
                -fq unmapped_sorted_1.fq \
                -fq2 unmapped_sorted_2.fq

    >>>

    output {
        File unmapped_aln_sam = "unmapped_aln.sam"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(unmapped_bam, "GB")*3 ) + " HDD"
        docker: docker
        cpu: cpus
        memory: "10GB"
    }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create the task BWA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
task BWA {
    input {
        
        File Virus_Reference
        File Human_Reference
        File GTF_Reference

        File unmapped_bam

        Int cpus
        Int preemptible
        String docker
        String sample_id
    }

    #~~~~~~~~~~~~~~~~~
    # Set the Prefix 
    #~~~~~~~~~~~~~~~~~
    String prefix = sample_id

    command <<<

        set -e

        #~~~~~~~~~~~~~~~
        # SORT
        #~~~~~~~~~~~~~~~
        samtools sort \
            -n ~{unmapped_bam} \
            -o unmapped_sorted.bam

        #~~~~~~~~~~~~~~~
        # Bedtools 
        #~~~~~~~~~~~~~~~

        bedtools bamtofastq -i . \
                unmapped_sorted.bam \
                -fq unmapped_sorted_1.fq \
                -fq2 unmapped_sorted_2.fq

    >>>

    output {
        File unmapped_aln_sam = "unmapped_aln.sam"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(unmapped_bam, "GB")*3 ) + " HDD"
        docker: docker
        cpu: cpus
        memory: "10GB"
    }
}






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

workflow VirTect {
    input {

        #~~~~~~~~~~~~
        # Sample ID
        #~~~~~~~~~~~~
        String sample_id
      
        #~~~~~~~~~~~~
        # FASTQ Files
        #~~~~~~~~~~~~
        File left
        File? right

        #~~~~~~~~~~~~
        # CPU count 
        #~~~~~~~~~~~~
        Int cpus = 1

        #~~~~~~~~~~~~
        # Directories 
        #~~~~~~~~~~~~
        File Virus_Reference
        File Human_Reference
        File GTF_Reference

        #~~~~~~~~~~~~
        # References
        #~~~~~~~~~~~~
        #File ref_fasta

        #~~~~~~~~~~~~
        # general runtime settings
        #~~~~~~~~~~~~
        Int preemptible = 2
        String docker = "brownmp/virtect:devel"

        

    }

    parameter_meta {
        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        cpus:{help:"CPU count"}
        docker:{help:"Docker image"}
    }

    call RunCutadapt{
        input:
            fastq1 = left,
            fastq2 = right,
            cpus = cpus, 

            preemptible = preemptible,
            docker = docker,
            sample_id = sample_id
    }

    call RunTopHat{
        input:
            fastq1          = RunCutadapt.fastq1_out,
            fastq2          = RunCutadapt.fastq2_out,

            Virus_Reference = Virus_Reference,
            Human_Reference = Human_Reference,
            GTF_Reference   = GTF_Reference,

            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id
    }

    call bam2fastq{
        input:
        
            unmapped_bam = RunTopHat.unmapped_bam,

            Virus_Reference = Virus_Reference,
            Human_Reference = Human_Reference,
            GTF_Reference   = GTF_Reference,

            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id
    }

    
}

