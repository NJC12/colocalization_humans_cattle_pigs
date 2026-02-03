#$ -S /bin/bash
#$ -e /net/home/nconnally/comparative_colocalization/fastenloc_simulations/eo
#$ -o /net/home/nconnally/comparative_colocalization/fastenloc_simulations/eo
#$ -l mf=3G
#$ -V
#$ -N dapg
#$ -pe serial 16

# On muhee compute cluster

eval "$(/net/apps/conda/miniconda3/bin/conda shell.bash hook)"
conda activate fact

# Arguments are:
# 1. The working directory /net/home/nconnally/comparative_colocalization/fastenloc_simulations/gwas_scaling_35.0_gtex_scaling_35.0_min_maf_0.01
# 2. Directory name of the category you're looking at
# 3. First line of phenotype file (*_pheno.sbams) you want to analyze (exclusive, starts at 0)
# 4. Last line of the above that you want to analyze
# 5. R2 cutoff (an argument used by dap-g)
# 6. First core to use. For some reason, on an SGE scheduler, dap-g grabs as many cpus as possible, so this limits it (in this example) to cpus 10-26.

# category should be cgtex, cgwas, hgtex, or hgwas
dir="$1"
category="$2"
startline="$3"
endline="$4"
ld_ctrl="$5"

# I have to assign specific CPUs because otherwise, on the cluster I'm using, dap-g grabs as many as it can regardless of the -t 4 flag
cpu_var="$6"
cpus=$(seq -s, $cpu_var $((cpu_var+15)))

cat_abbv=$(echo $category | sed -r 's/_.*//g')
cd $dir/$cat_abbv
mkdir -p tmp_sbams outputs logs

for ph in $(awk -v start=$startline -v end=$endline 'NR > start && NR <= end {print $2}' ${category}_pheno.sbams); do
    if [ -f "outputs/${category}_$ph.dapg.out" ]; then
        echo "File for $tr exists--skipping"
        continue
    fi
    position=$(echo $ph | sed 's/tr//g')
    echo $ph -- position $position
    awk -v trait=$ph '$2 == trait' ${category}_pheno.sbams > tmp_sbams/${category}_$ph.sbams
    # cat ${category}_geno.sbams >> tmp_sbams/${category}_$ph.sbams
    awk -v pos=$position 'BEGIN{OFS="\t"} {gsub(/snp/, "", $2)} $2+0 >= pos-1e6 && $2+0 <= pos+1e6 {print $0}' ${category}_geno.sbams >> tmp_sbams/${category}_$ph.sbams
    cat ../pca/$cat_abbv.pca >> tmp_sbams/${category}_$ph.sbams
    taskset -c $cpus dap-g -d tmp_sbams/${category}_$ph.sbams --output_all -t 16 -msize 5 -o outputs/${category}_$ph.dapg.out -l logs/${category}_$ph.dapg.out -ld_control ld_ctrl

    rm tmp_sbams/${category}_$ph.sbams
done