#!/usr/bin/env nextflow


/*
 * Pipeline Metagenomics, conda ambient Metagenomics_Nextflow The conda ambient incled the following packages::: The links could be broken
    * Bbmap https://anaconda.org/bioconda/bbmap
    * Bowtie2 https://anaconda.org/bioconda/bowtie2
    * Samtools https://anaconda.org/bioconda/samtools
    * Diammond https://anaconda.org/bioconda/diamond
    * Famli https://github.com/FredHutch/FAMLI Not available on conda
    * Upimapi https://anaconda.org/bioconda/kraken2
 */


 nextflow.enable.dsl=2 // Nextflow version

params.fq1 = "$HOME/../*1.fq.gz" // Can be modified in the script with the --fq1 <value> command | Put read1
params.fq2 = "$HOME/../*2.fq.gz" // Can be modified in the script with the --fq2 <value> command | Put read2
params.db_human = "./GRCh38_latest_genomic.fna"// Can be modified in the script with the --human <value> command | Put bowtie index human genome
params.db_uniref = "./Kraken_raw/"// Can be modified in the script with the --db_gut <value> command | The Uniref 90 uniprot database
params.otuput_dir = "$HOME/../" // Can be modified in the script with the --otuput_dir <value> command | Put OTU table dir output


/*
 * Check the quality of raw sequences
 */

 process BBduk {

    input:
        path read1
        path read2
    output:
        file("*_1.fq.tr.gz")
        file("*_2.fq.tr.gz")
        file("*.bbduk.log") 

     script:
        def prefix = read1.toString().replaceFirst("_1\\.fq\\.gz\$", "")
        def raw = "in1=${read1} in2=${read2}"
        def trimmed  = "out1=${prefix}_1.fq.tr.gz out2=${prefix}_2.fq.tr.gz"

            """
                bbduk.sh \\
                $raw \\
                $trimmed \\
                qtrim=r trimq=10 minlen=100 \\
                &> ${prefix}.bbduk.log
            """

 }

/*
* Remove the human DNA sequences
*/

 process DNA_clean {

     input:
        path file1
        path file2
        val DB 
    output:
        file "*.R1.fq.gz"
        file "*.R2.fq.gz"
        file "*.Bowtie.log"
    script:
        def prefix = file1.toString().replaceFirst("_1\\.fq\\.tr\\.gz\$", "")
        def entry_bowtie = "-x ${DB} -1 ${file1} -2 ${file2}"
        """
            # Math with the human chromosome
                bowtie2 $entry_bowtie | samtools  view -Sb - > alignment.bam 2 > ${prefix}.Bowtie.log
            # Flag all the mathed Human_DNA
                samtools view -b -f 12 alignment.bam > archivo.filt.bam
            # Create both reads 
                samtools fastq archivo.filt.bam -1 ${prefix}.R1.fq.gz -2 ${prefix}.R2.fq.gz
        """
 }

/*
* Assing the protein ID in especified database
*/

 process Protein_assingment {
     input:
        path file1
    publishDir  "${params.otuput_dir}", mode: 'copy' // Check if it not needed
    output:
        file "*.Uniprot.tsv"
    script:
        def prefix = file1.toString().replaceFirst(".R1\\.fq\\.gz\$", "")
        def entry_diamond = "${file1} ${file2}"
        """
            # Math with the human chromosome
                cat ${entry_diamond} |
                diamond blastx --out $name.diamond.gz --compress 1 
                  --threads 15 --db ${DB} 
                  --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen --query-cover 90 --id 80 --top 0 --query-gencode 11 --unal 0
            # Clean the otuput
                 famli filter 
                   --threads 15 --qseqid-ix 0 --sseqid-ix 1 --sstart-ix 6 --send-ix 7 --bitscore-ix 9 --slen-ix 11
            # Create the table
              jq -r '.[] | .id, .nreads' | paste - - 
        """
 }
/*
* Annotate the proteins
*/
 process Protein_anotation {
     input:
        path file1
    publishDir  "${params.otuput_dir}", mode: 'copy' // Check if it not needed
    output:
        file "diamond.upimapi"
        """
            # Annotate the preotins
            cut -f1 *${file1}.tsv | sort -u
            upimapi.py -o diamond.upimapi
        """
 }




/*
* Workflow (Complete proces)
*/
workflow {
      sequences = [file(params.fq1), file(params.fq2)]
    println "Performing Quality control and trimming from $sequences" 

    // Run the BBduk process
    Output_bbduk = BBduk(sequences[0], sequences[1])
    
    // Check the human database and create the database variable
    println "The Human Database selected is $params.db_uniref" 
    Databases = [file(params.db_human), file(params.db_uniref)]
    
    // The outputs of BBduk are automatically connected to DNA_clean
    Output_bowtie = DNA_clean(Output_bbduk[0], Output_bbduk[1], Databases[0])

    // Obtain proteins IDs for each read
    Protein_assingment(Output_bowtie[0], Output_bowtie[1], Databases[1])
    
    // Create the annotation for each ID
    println "The output folder selected is $params.otuput_dir"
    Protein_anotation(Output_Protein_assingment[0])
}
