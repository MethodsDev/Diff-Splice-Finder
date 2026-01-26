version 1.0

workflow CountIntrons {
    input {
        File bam_file
        File bam_index
        File genome_fasta
        String docker = "methodsdev/diff-splice-finder:latest"
        Int cpu = 4
        Int memory_gb = 8
        Int disk_gb = 100
    }

    call CountIntronsFromBam {
        input:
            bam_file = bam_file,
            bam_index = bam_index,
            genome_fasta = genome_fasta,
            docker = docker,
            cpu = cpu,
            memory_gb = memory_gb,
            disk_gb = disk_gb
    }

    output {
        File intron_counts = CountIntronsFromBam.intron_counts
    }
}

task CountIntronsFromBam {
    input {
        File bam_file
        File bam_index
        File genome_fasta
        String docker
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    String output_filename = basename(bam_file, ".bam") + ".intron_counts.tsv"

    command <<<
        set -e
        
        # Run the intron counting script
        count_introns_from_bam.py \
            --bam ~{bam_file} \
            --genome_fa ~{genome_fasta} \
            > ~{output_filename}
    >>>

    output {
        File intron_counts = output_filename
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
    }
}
