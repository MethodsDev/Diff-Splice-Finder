version 1.0

workflow CountIntrons {
    input {
        String sample_id
        File bam_file
        File bam_index
        File genome_fasta
        File? target_introns
        Int site_depth_window_radius = 10
        Int min_mapping_quality = 60
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
            target_introns = target_introns,
            site_depth_window_radius = site_depth_window_radius,
            min_mapping_quality = min_mapping_quality,
            docker = docker,
            cpu = cpu,
            memory_gb = memory_gb,
            disk_gb = disk_gb
    }

    output {
        File introns = CountIntronsFromBam.introns
        File intron_counts = CountIntronsFromBam.intron_counts
    }
}

task CountIntronsFromBam {
    input {
        String sample_id
        File bam_file
        File bam_index
        File genome_fasta
        File? target_introns
        Int site_depth_window_radius
        Int min_mapping_quality
        String docker
        Int cpu
        Int memory_gb
        Int disk_gb
    }

    String output_filename = "~{sample_id}.introns.gz"

    command <<<
        set -euo pipefail
        
        # Run the intron counting script and gzip output
        count_introns_from_bam.py \
            --bam ~{bam_file} \
            --genome_fa ~{genome_fasta} \
            ~{if defined(target_introns) then "--target_introns " + select_first([target_introns]) else ""} \
            --site_depth_window_radius ~{site_depth_window_radius} \
            --min_mapping_quality ~{min_mapping_quality} \
            | gzip -c > ~{output_filename}
    >>>

    output {
        File introns = output_filename
        File intron_counts = output_filename
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
    }
}
