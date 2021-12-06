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
task RunVirTect {
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

        ls -lt 

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

        cp ~{GTF_Reference} .

        #~~~~~~~~~~~~~~~
        # VirTect
        #~~~~~~~~~~~~~~~
        python3 /usr/local/src/VirTect/VirTect.py \
                -t 12 \
                -1 ~{fastq1} \
                -2 ~{fastq2} \
                -o ~{prefix} \
                -ucsc_gene ~{GTF_Reference} \
                -index GRCh38.genome \
                -index_vir ~{Virus_Reference} \
                -d 200 \
                --n_thread ~{cpus}


        samtools depth ~{prefix}/unmapped_aln_sorted.bam | awk '{if ($3>=5) print $0}'| awk '{ if ($2!=(ploc+1)) {if (ploc!=0){printf( "%s %d-%d\n",$1,s,ploc);}s=$2} ploc=$2; }' \
            > ~{prefix}/continuous_region.txt


        python <<CODE

            print ("The continous length")
            file =open(out+"/continuous_region.txt", "r")
            out_put =open(out+"/Final_continous_region.txt", "w")
            
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
                

            final_output=open(out+"/Final_continous_region.txt",'r')
            if (os.fstat(final_output.fileno()).st_size) >0:
                print ("----------------------------------------Note: The sample may have some real virus :(-----------------------------------------------------")
                headers = 'virus transcript_start transcript_end'.split()
                for line in fileinput.input([out+'/Final_continous_region.txt'], inplace=True):
                    if fileinput.isfirstline():
                        print '\t'.join(headers)
                    print line.strip()
            else:
                print ("----------------------------------------Note: There is no real virus in the sample :)-----------------------------------------------------")

        CODE


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

