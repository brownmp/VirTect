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


        python /usr/local/src/VirTect/VerTect_cutadapt.py \
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
task RunVirTect {
    input {
        File fastq1
        File? fastq2
        
        File Virus_Reference
        #File Human_Reference
        File GTF_Reference

        ref_fasta = ref_fasta

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

        ls -lt 

        # special case for tar of fastq files
        if [[ "~{fastq1}" == *.tar.gz ]] ; then
            mkdir fastq
            tar -I pigz -xvf ~{fastq1} -C fastq
            fastqs=$(find fastq -type f)
            fastq1=$fastqs[0]
            fastq2=$fastqs[1]
        fi


        #~~~~~~~~~~~~~~~
        # VirTect
        #~~~~~~~~~~~~~~~
        python3 /usr/local/src/VirTect/VirTect.py \
                -t 12 \
                -1 ~{fastq1} \
                -2 ~{fastq2} \
                -o ~{prefix} \
                -ucsc_gene ~{GTF_Reference} \
                -index ~{ref_fasta} \
                -index_vir ~{Virus_Reference} \
                -d 200 \
                --n_thread ~{cpus}


        samtools depth ~{prefix}/unmapped_aln_sorted.bam | awk '{if ($3>=5) print $0}'| awk '{ if ($2!=(ploc+1)) {if (ploc!=0){printf( "%s %d-%d\n",$1,s,ploc);}s=$2} ploc=$2; }' \
            > ~{prefix}/continuous_region.txt

    >>>

    output {

    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fastq1, "GB")*3 ) + " HDD"
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
        #File Human_Reference
        File GTF_Reference

        #~~~~~~~~~~~~
        # References
        #~~~~~~~~~~~~
        File ref_fasta

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

    call RunVirTect{
        input:
            fastq1          = RunCutadapt.fastq1_out,
            fastq2          = RunCutadapt.fastq2_out,
            Virus_Reference = Virus_Reference,
            #Human_Reference = Human_Reference,
            GTF_Reference   = GTF_Reference,

            ref_fasta = ref_fasta
            

            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id
    }

    
}

