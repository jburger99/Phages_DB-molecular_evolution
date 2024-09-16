process get_fastas {
    conda "conda-forge::r-tidyverse=2.0.0"
    conda "conda-forge::r-biocmanager=1.30.24"
    conda "conda-forge::r-fuzzyjoin=0.1.6"
    conda "conda-forge::r-rlist=0.4.6.2"
    input:
    val x
    output:
    path "${x}*.fas"
    script:
    """
    #!/usr/bin/env Rscript
packages = c("tidyverse", "BiocManager", "fuzzyjoin")

## load or install & load all packages
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

packages <- c("Biostrings", "seqinr")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


## import phagesDB fasta file

url <- "https://phagesdb.org/media/Actinobacteriophages-All.fasta"

suppressWarnings(all_phages <- readDNAStringSet(url))

## Create a dataframe from fasta file

all_phages_names <- names(all_phages)
all_phages_sequences <- paste(all_phages)
all_phages_df <- data.frame(all_phages_names, all_phages_sequences)

## import phage names and clusters from phagesdb

names_clusters <- read_tsv("https://phagesdb.org/data/?set=seq&type=simple")

all_phages_df_final <- names_clusters %>%
  fuzzyjoin::regex_left_join(all_phages_df, ., by = c(all_phages_names = "Phage Name")) %>% 
  filter(!is.na(`Phage Name`)) %>% 
  reframe(name = `Phage Name`, Cluster, Subcluster, seq = all_phages_sequences)%>% 
  filter(!duplicated(name)) 

all_phages_df_final[all_phages_df_final == "None"] <- NA

all_phages_df_final <- all_phages_df_final %>% 
  mutate(Subcluster = coalesce(Subcluster, Cluster))

## get_phams script

phams = ${x}

phams = as_tibble(phams) |>
  reframe(Pham = as.numeric(value),
          URL = paste("https://phagesdb.org/phams/genelist/",
                      as.numeric(value),
                      "/", 
                      sep = ""))



all_phages <- all_phages_df_final

## Read in the Gene Data acquired from PhagesDB
phages_targets <- data.frame() 

for (i in phams[,2]){
  pham <- read_tsv(i)
  phages_targets <- rbind(phages_targets, pham)
}

colnames(phages_targets)[1] <- "name"

phages_targets <- phages_targets |> 
  reframe(name, `Locus Tag`, Start, Stop, pham)

## Join Gene Data with Phage Data

phages_target_joined <- inner_join(all_phages, 
                                   phages_targets, by = "name") 
phages_target_joined <- phages_target_joined |> 
  mutate(count = Stop - Start) |>
  filter(`Locus Tag` != "None")
  

## get_seqs_by_cluster Script

clusters <- phages_target_joined |> 
  group_by(Subcluster) |> 
  summarise(n = n()) %>% 
  filter(n >= 3)



## For loop to extra the desired gene sequence from the Phage Seq Data


for (i in 1:nrow(clusters)){
  
  phages_target_joined_loop <- phages_target_joined |> 
    filter(Subcluster == as.character(clusters[i,1]))
  
  seqs <- list(1:nrow(phages_target_joined_loop))
  
  
  
  seqs_target <- list()
  
  for(j in seqs){
    value <- as.numeric(j)
    start <-  as.numeric(phages_target_joined_loop[value,6] + 3)
    end <- as.numeric(phages_target_joined_loop[value,7] - 3)
    
    
    gene_seq <- str_sub(phages_target_joined_loop[j,4], start = start, end = end)
    
    
    seqs_target <- append(seqs_target, gene_seq)
  
  
  }
  
 
  ## Creation of FASTA File
  seq_table <- data.frame(phages_target_joined_loop[,1], as.character(seqs_target)) %>% 
    reframe(names = `phages_target_joined_loop...1.`, 
            seqs = lapply(`seqs_target`, DNAString)) %>% 
    filter(seqs_target != "")
  names <- seq_table["names"]
  write.fasta(sequences = as.list(seq_table[,2]), names = names[,1],
              file.out = paste(${x}, "_", clusters[i,1], "_Seqs.fas", sep = ""))
  
}
  

    """
}

process transeq {
  conda "bioconda::emboss=6.6.0"

  input:

  path x

  output:

  path '*_amino.fas'

  script:
  """
  #!/usr/bin/env python
import subprocess

file_type = ".fas"
string = "${x}"
files = list(string.split())
for filename in files:
    if filename.endswith(file_type):
        amino_file = filename.rstrip(".fas") + "_amino.fas"
        transeq = "transeq -sequence " + filename + " -outseq " + amino_file
        subprocess.run(transeq, shell=True)

"""


}

process clustalo_align {
  conda "bioconda::clustalo=1.2.4"
  input:
  path x

  output:
  path '*_aligned.fas'

  script:
"""
#!/usr/bin/env python
import subprocess

file_type = ".fas"
string = "${x}"
files = list(string.split())
for filename in files:
    if filename.endswith(file_type):
        aligned_file = filename.rstrip("amino_aligned.fas") + "_aligned.fas"
        clustal = "clustalo -i " + filename + " -o " + aligned_file
        subprocess.run(clustal, shell=True)

"""

}

process tranalign {

  conda "bioconda::clustalo=1.2.4"

  input:

  path x
  path y

  output:

  path '*ntaligned.fas'

  script:
  """
  #!/usr/bin/env python
  import subprocess

  file_type = "aligned.fas"
  string = "${x}"
  string2 = "${y}"
  files = list(string.split())
  for filename in files:
      if filename.endswith(file_type):
        nt_file = filename.rstrip("amino_alinged.fas") + "s.fas"
        command = "tranalign -asequence " + nt_file + " -bsequence " + filename + " -outseq " + filename + "ntaligned.fas"

        subprocess.run(command, shell=True)

  """

}

process modeltest {
  conda "bioconda::modeltest-ng=0.1.7"

  input:
  path x

  output: 
  path '*.out'

  script:
  """
  #!/usr/bin/env python
  import subprocess

  file_type = "aligned.fas"
  string = "${x}"
  files = list(string.split())
  for filename in files:
      if filename.endswith(file_type):
        command = "modeltest-ng -i " + filename + " -o " + filename + ".out"

        subprocess.run(command, shell=True)

  """
}

process iqtree {
  input:
  path x
  path y

  output:
  path '*'

  script:
  """
  #!/usr/bin/env python

import re
import subprocess

fastas = "${y}"

def find_command(filename):

    with open(filename, "r") as file:
     lines = file.readlines()
     results = list()
     for row in lines:
        command = re.compile("iqtree -s .* -m .*")
        match = command.search(row)
        if match:
            results.append(match)
     clean_command = re.sub(".*match='", "", str(results[2]))
     final_command = clean_command.rstrip("'>") 

     subprocess.run(final_command, shell=True)
     

def run_iqtree(path):
   file_type = ".out.out"
   for file in path:
      if file.endswith(file_type):
            find_command(file)
         
      
string = "${x}"
path = list(string.split()) 

run_iqtree(path)

"""

}

process meme {
  conda "bioconda::hyphy=2.5.6"

  input:
  path x
  path y

  output:

  path '*_MEME.json'

  script:
  """
  #!/usr/bin/env python

  import subprocess

  treefiles = "${y}"

  file_type = "aligned.fas"
  string = "${x}"
  files = list(string.split())
  for filename in files: 
      if filename.endswith(file_type):
          command = "hyphy meme --alignment " + filename + " --tree " + filename + ".treefile" + " --code Universal --output " + filename + "_MEME.json"

          subprocess.run(command, shell=True)

  """
  
}

workflow{
    def value = params.pham
    get_fastas(value)
    transeq(get_fastas.out)
    clustalo_align(transeq.out)
    tranalign(clustalo_align.out, get_fastas.out)
    modeltest(tranalign.out)
    iqtree(modeltest.out, tranalign.out)
    meme(tranalign.out, iqtree.out)

}