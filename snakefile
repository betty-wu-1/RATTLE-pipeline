rule pc:
	input: "replicates/{sample}.fq"
	output: "pc_output/1.{sample}.fq"
	shell: "~/softwares/Porechop/porechop-runner.py -i {input} -o {output}"

rule minimap2:
        input: "pc_output/1.{sample}.fq"
        output: "minimap2_output/1.2.{sample}.sam"
        shell: "~/minimap2/minimap2 -f -a -ax splice --junc-bonus 1 -k14 --secondary=no --junc-bed -uf -k14 -t24 /g/data/xc17/betty/human_reference_genome/hg38.fa {input}> {output}"

rule sorted_bam:
	input: "minimap2_output/1.2.{sample}.sam"
        output: "minimap2_output/1.2.{sample}.sorted.bam"
	shell: """
		module load samtools
		samtools view -bS {input} | samtools sort -@ 48 > {output}
		"""
	
rule unmapped_sam:
	input: "minimap2_output/1.2.{sample}.sorted.bam"
	output: "minimap2_output/unmapped.{sample}.sam"
	shell: """
                module load samtools
		samtools view -f 4 {input}> {output}
		"""

rule sam_fastq:
	input: "minimap2_output/unmapped.{sample}.sam"	
	output: "minimap2_output/1.2.{sample}.fq"
	shell: r"""cat {input} | grep -v ^@ | awk '{{print "@"$1"\n"$10"\n+\n"$11}}' > {output}"""

rule nanofilt:
	input: "minimap2_output/1.2.{sample}.fq"
	output: "nanofilt_output/1.2.3.{sample}.fq"
	shell: """cat {input}| NanoFilt -q 7 > {output}
		mkdir -p rattle_output/{wildcards.sample}_results"""
 
rule rattle_cluster:
	input: "nanofilt_output/1.2.3.{sample}.fq"
	output: "rattle_output/{sample}_results/clusters.out"
	shell: """~/softwares/RATTLE/rattle cluster -i {input} -t 24 -o rattle_output/{wildcards.sample}_results -B 0.5 -b 0.3 -f 0.2 --fastq --rna"""

rule rattle_correct:
	input: 
		clusters="rattle_output/{sample}_results/clusters.out",
		fasta="nanofilt_output/1.2.3.{sample}.fq"
	output: "rattle_output/{sample}_results/consensi.fq"
	shell: "~/softwares/RATTLE/rattle correct -i {input.fasta} -c {input.clusters} -o rattle_output/{wildcards.sample}_results -t 24"

rule rattle_polish:
	input: "rattle_output/{sample}_results/consensi.fq"
	output: "rattle_output/{sample}_results/transcriptome.fq"
	shell: "~/softwares/RATTLE/rattle polish -i {input} -o rattle_output/{wildcards.sample}_results -t24 --rna"
