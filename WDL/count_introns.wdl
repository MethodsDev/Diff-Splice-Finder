version 1.0

workflow CountIntrons {
    input {
        String sample_id
        File bam_file
        File bam_index
        File genome_fasta
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/diff-splice-finder/diff-splice-finder"
        Int cpu = 4
        Int memory_gb = 8
        Int disk_gb = 100
    }

    call CountIntronsFromBam {
        input:
            sample_id = sample_id,
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
        String sample_id
        File bam_file
        File bam_index
        File genome_fasta
        String docker
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    String output_filename = "~{sample_id}.introns"

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
