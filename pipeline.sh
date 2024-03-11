# Define a variable where is cloned the github repository
gitHubDirectory=/home/ldelisle/softwares/scriptsAroundBaredSCForDarbellayEtAl2024

# On the cluster, define a working directory
mkdir -p /scratch/ldelisle/baredSC_Darbellay
cd /scratch/ldelisle/baredSC_Darbellay

# Create a conda environment
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda create -n baredSC baredsc=1.1.2

# Split the full table by merged.cluster
inputTable=${gitHubDirectory}/inputs/sc_metadata.txt
cat ${inputTable} | awk -F "\t" -v output="merged.cluster" -v OFS="\t" '
NR == 1{
  for(k = 1; k <= NF; k++){
    if($k == "merged.clusters"){
      i = k
    }
  }
  header=$0
}
NR>1{
  cluster = $i
  if(! (cluster in seen)){
    # Print the header
    print header > output"_"cluster".txt"
    seen[cluster] = 1
  }
  print $0 > output"_"cluster".txt"
}'

# Split the mesenchyme table by cluster
inputTable=${gitHubDirectory}/inputs/sc_mesenchyme_metadata.txt
cat ${inputTable} | awk -F "\t" -v output="cluster" -v OFS="\t" '
NR == 1{
  for(k = 1; k <= NF; k++){
    if($k == "original_cluster_identity"){
      i = k
    }
  }
  header=$0
}
NR>1{
  cluster = $i
  if(! (cluster in seen)){
    # Print the header
    print header > output"_"cluster".txt"
    seen[cluster] = 1
  }
  print $0 > output"_"cluster".txt"
}'

# Prepare table for 2D
pathForTable="${gitHubDirectory}/tables/table_baredSC_2d.txt"
mkdir -p $(dirname $pathForTable)
if [ -e $pathForTable ]; then
    rm $pathForTable
fi
for f in cluster_*.txt; do
    echo -e "${f}\tCol2a1\teGFP-SV40pA\t6\t4" >> $pathForTable
done
# I launch the parallel baredSC in 2D for mesenchyme clusters:
sbatch --array 1-8 --chdir $PWD/ ${gitHubDirectory}/sbatch_scripts/sbatch_baredSC_2d.sh ${pathForTable} $PWD/baredSC_2d/

# Then we generate the table for parallel baredSC in 1d:
pathForTable="${gitHubDirectory}/tables/table_baredSC_1d.txt"
mkdir -p $(dirname $pathForTable)
if [ -e $pathForTable ]; then
  rm $pathForTable
fi
nsim=0
for f in cluster_*.txt; do
  for gene in Col2a1 eGFP-SV40pA ; do
    echo -e "${f}\t${gene}\t6" >> $pathForTable
    nsim=$((nsim + 1))
  done
done
for f in merged.cluster_*.txt; do
  for gene in Col2a1 eGFP-SV40pA ; do
    echo -e "${f}\t${gene}\t6" >> $pathForTable
    nsim=$((nsim + 1))
  done
done
sbatch --array 1-${nsim} --chdir $PWD/ ${gitHubDirectory}/sbatch_scripts/sbatch_baredSC_1d.sh ${pathForTable} $PWD/baredSC_1d/

# Once everything is done:
module load gcc/11.3.0 r/4.1.3
Rscript ${gitHubDirectory}/r_scripts/MCMC_1d_plots_per_merged.cluster.R $PWD/ "baredSC_1d" ${gitHubDirectory}/tables/table_baredSC_1d.txt ${gitHubDirectory}/plots/merged.clusters
Rscript ${gitHubDirectory}/r_scripts/MCMC_1d_plots_per_cluster.R $PWD/ "baredSC_1d" ${gitHubDirectory}/tables/table_baredSC_1d.txt ${gitHubDirectory}/plots/clusters
Rscript ${gitHubDirectory}/r_scripts/MCMC_2d_plots.R $PWD/ "baredSC_2d" ${gitHubDirectory}/tables/table_baredSC_2d.txt ${gitHubDirectory}/plots/2d
