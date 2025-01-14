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
        if [[ "~{fastq1}" == *.tar.gz ]]; then
            mkdir fastq
            tar -I pigz -xvf ~{fastq1} -C fastq
            fastqs=$(find fastq -type f)
            fastq1=$fastqs[0]
            fastq2=$fastqs[1]
        fi


        python3 /usr/local/src/VirTect/VerTect_cutadapt.py \
                -1 $fastq1 \
                -2 $fastq2 \
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
        disks: "local-disk " + ceil(size(fastq1, "GB") * 4.5 ) + " HDD"
        docker: docker
        cpu: 1
        memory: "30GB"
    }
}




####################################
# create the task RunVirTect
####################################
task RunTopHat {
    input {
        File fastq1
        File? fastq2
        
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

        # Untar the references  
        tar -xvf ~{Human_Reference}

        # special case for tar of fastq files
        if [[ "~{fastq1}" == *.tar.gz ]] ; then
            mkdir fastq
            tar -I pigz -xvf ~{fastq1} -C fastq
            fastqs=$(find fastq -type f)
            
            #~~~~~~~~~~~~~~~
            # TopHat
            #~~~~~~~~~~~~~~~
            tophat2 \
                -o ~{prefix} \
                -p ~{cpus} \
                -G ~{GTF_Reference} \
                GRCh38.genome \
                $fastqs

        else
            #~~~~~~~~~~~~~~~
            # TopHat
            #~~~~~~~~~~~~~~~
            tophat2 \
                -o ~{prefix} \
                -p ~{cpus} \
                -G ~{GTF_Reference} \
                GRCh38.genome \
                ~{fastq1} ~{fastq2}

        fi
    >>>

    output {
        File unmapped_bam = "~{prefix}/unmapped.bam"
        File accepted_hits_bam = "~{prefix}/accepted_hits.bam"
        File junctions_bed = "~{prefix}/junctions.bed"
        File deletions_bed = "~{prefix}/deletions.bed"
        File insertions_bed = "~{prefix}/insertions.bed"
        File align_summary_txt = "~{prefix}/align_summary.txt"
        File prep_reads_info = "~{prefix}/prep_reads.info"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fastq1, "GB")*6 + size(GTF_Reference, "GB") + size(Human_Reference, "GB")*3 + 50) + " HDD"
        docker: docker
        cpu: cpus
        memory: "150GB"
    }
}



####################################
# create the task bam2fastq
####################################

task bam2fastq {
    input {
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
            -o ~{prefix}_unmapped_sorted.bam

        #~~~~~~~~~~~~~~~
        # Bedtools 
        #~~~~~~~~~~~~~~~
        bedtools bamtofastq \
            -i ~{prefix}_unmapped_sorted.bam \
            -fq ~{prefix}_unmapped_sorted_1.fq \
            -fq2 ~{prefix}_unmapped_sorted_2.fq

    >>>

    output {
        File unmapped_sorted_bam = "~{prefix}_unmapped_sorted.bam"
        File unmapped_sorted_1 = "~{prefix}_unmapped_sorted_1.fq"
        File unmapped_sorted_2 = "~{prefix}_unmapped_sorted_2.fq"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(unmapped_bam, "GB")*6 + 50) + " HDD"
        docker: docker
        cpu: 1
        memory: "100GB"
    }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create the task BWA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
task BWA {
    input {
        
        File unmapped_sorted_1
        File unmapped_sorted_2

        File Virus_Reference
        File GTF_Reference

        Int cpus
        Int preemptible
        String docker
        String sample_id
        String Virus_Fasta_Prefix
    }

    #~~~~~~~~~~~~~~~~~
    # Set the Prefix 
    #~~~~~~~~~~~~~~~~~
    String prefix = sample_id

    command <<<

        set -e

        # Untar the references  
        tar -xvf ~{Virus_Reference}


        #~~~~~~~~~~~~~~~
        # BWA
        #~~~~~~~~~~~~~~~

        bwa mem \
            ~{Virus_Fasta_Prefix}.fasta \
            ~{unmapped_sorted_1} \
            ~{unmapped_sorted_2} \
            > ~{prefix}_unmapped_aln.sam

    >>>

    output {
        File unmapped_aln_sam = "~{prefix}_unmapped_aln.sam"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(unmapped_sorted_1, "GB")*3 + size(GTF_Reference, "GB") + size(Virus_Reference, "GB")*3 ) + " HDD"
        docker: docker
        cpu: cpus
        memory: "10GB"
    }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create the task Virus Detection
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
task VirusDetection {
    input {
        
        File unmapped_aln_sam

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
        # Samtools
        #~~~~~~~~~~~~~~~
        samtools view -Sb \
                -h ~{unmapped_aln_sam} \
                > ~{prefix}_unmapped_aln.bam


        samtools view \
            ~{prefix}_unmapped_aln.bam | cut -f3 | sort | uniq -c | awk '{if ($1>=400) print $0}' \
            > ~{prefix}_unmapped_viruses_count.txt

    >>>

    output {
        File unmapped_aln_bam = "~{prefix}_unmapped_aln.bam"
        File unmapped_viruses_count = "~{prefix}_unmapped_viruses_count.txt"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(unmapped_aln_sam, "GB")*4 ) + " HDD"
        docker: docker
        cpu: cpus
        memory: "100GB"
    }
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ContinuousRegion
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
task ContinuousRegion {
    input {
        
        File unmapped_aln_bam

        Int Distance

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
        # Samtools
        #~~~~~~~~~~~~~~~
        samtools sort \
            ~{unmapped_aln_bam} \
            -o ~{prefix}_unmapped_aln_sorted.bam

        samtools depth ~{prefix}_unmapped_aln_sorted.bam | awk '{if ($3>=5) print $0}'| awk '{ if ($2!=(ploc+1)) {if (ploc!=0){printf( "%s %d-%d\n",$1,s,ploc);}s=$2} ploc=$2; }' \
            > ~{prefix}_continuous_region.txt


        #~~~~~~~~~~~~~~~
        # Python 
        #~~~~~~~~~~~~~~~
        python <<CODE

        import sys
        import argparse
        import os
        import subprocess
        import os.path
        import fileinput

        distance = ~{Distance}

        print("The continous length")
        file =open("~{prefix}_continuous_region.txt", "r")
        out_put =open("~{prefix}_Final_continous_region.txt", "w")
        
        if (os.fstat(file.fileno()).st_size) >0:
                for i in file.readlines():
                    i1=i.split()[0]
                    i2=i.split()[1]
                    j1=i2.split("-")
                    j2=int(j1[1])-int(j1[0])

                    if j2 >= distance:
                        j3=i1 + "\t" +  str(j1[0]) + '\t' +  str(j1[1])
                        out_put.write('%s\n' % j3)
                    else:
                        pass
        else:
            pass 
        out_put.close()
            

        final_output=open("~{prefix}_Final_continous_region.txt",'r')
        if (os.fstat(final_output.fileno()).st_size) >0 :
            print("----------------------------------------Note: The sample may have some real virus :(-----------------------------------------------------")
            headers = 'virus transcript_start transcript_end'.split()
            for line in fileinput.input(['~{prefix}_Final_continous_region.txt'], inplace=True):
                if fileinput.isfirstline():
                    print('\t'.join(headers))
                print(line.strip())
        else:
            print("----------------------------------------Note: There is no real virus in the sample :)-----------------------------------------------------")

        CODE

    >>>

    output {
        File unmapped_aln_sorted_bam = "~{prefix}_unmapped_aln_sorted.bam"
        File continuous_region = "~{prefix}_continuous_region.txt"
        File Final_continuous_region = "~{prefix}_Final_continous_region.txt"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(unmapped_aln_bam, "GB")*2 ) + " HDD"
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

        # Prefix for the virus fasta file
        String Virus_Fasta_Prefix
      
        #~~~~~~~~~~~~
        # FASTQ Files
        #~~~~~~~~~~~~
        File left
        File? right

        #~~~~~~~~~~~~
        # CPU count 
        #~~~~~~~~~~~~
        Int cpus = 10

        #~~~~~~~~~~~~
        # Directories 
        #~~~~~~~~~~~~
        File Virus_Reference
        File Human_Reference
        File GTF_Reference

        # Distacne value 
        Int Distance = 200

        #~~~~~~~~~~~~
        # Cutting the adaptors 
        #~~~~~~~~~~~~
        #Boolean cut_adaptors = false

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

    ##~~~~~~~~~~~~
    ## Cut the adaptors 
    ##~~~~~~~~~~~~
    #call RunCutadapt{
    #    input:
    #        fastq1 = left,
    #        fastq2 = right,
    #        cpus = cpus, 
    #
    #        preemptible = preemptible,
    #        docker = docker,
    #        sample_id = sample_id
    #}

    call RunTopHat{
        input:
            fastq1          = left,
            fastq2          = right,

            #Virus_Reference = Virus_Reference,
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

            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id
    }

    call BWA{
        input:
        
            unmapped_sorted_1 = bam2fastq.unmapped_sorted_1,
            unmapped_sorted_2 = bam2fastq.unmapped_sorted_2,

            Virus_Reference = Virus_Reference,
            #Human_Reference = Human_Reference,
            GTF_Reference   = GTF_Reference,

            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id,
            Virus_Fasta_Prefix = Virus_Fasta_Prefix
    }

    call VirusDetection {
        input:
            unmapped_aln_sam = BWA.unmapped_aln_sam,

            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id
    }

    call ContinuousRegion {
        input:
            unmapped_aln_bam = VirusDetection.unmapped_aln_bam,

            Distance = Distance,

            cpus            = cpus,
            preemptible     = preemptible,
            docker          = docker,
            sample_id       = sample_id
    }

    
}

