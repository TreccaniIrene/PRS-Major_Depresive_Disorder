#####################################
#                                   #
#            BASE                   #
#                                   #
#####################################

# Download LifetimeMDD database
wget https://figshare.com/ndownloader/files/21356610



#####################################
#                                   #
#           TARGET DATA QC          #
#                                   #
#####################################

# PLINK quality control (QC) steps on the "Target" dataset
plink --bfile Target \
      --snps-only just-acgt \
      --biallelic-only \
      --maf 0.01 \    # Removes all SNPs with minor allele frequency 
      --hwe 1e-6 \    # Removes SNPs with low P-value from the Hardy-Weinberg Equilibrium Fisher's exact or chi-squared test
      --geno 0.1 \    # Excludes SNPs that are missing in a high fraction of subjects
      --mind 0.1 \    # Excludes individuals who have a high rate of genotype missingness
      --allow-no-sex \
      --make-bed \
      --out Target_QC
      
# Performing LD pruning
# Using a window size of 50 variants, sliding across the genome with a step size of 5 variants at a time
# Filtering out any SNPs with a linkage disequilibrium (LD) r2 higher than 0.2
plink --bfile Target_QC \
      --indep-pairwise 50 5 0.2 \
      --out Target_QC
      
# Filtering out SNPs based on LD pruning results, calculating heterozygosity rates and saving the results as "pgc_clear"
plink --bfile Target_QC \
      --extract Target_QC.prune.in \   # Using the list of SNPs obtained from LD pruning
      --keep pgc_QC.fam \           # Keeping the corresponding individuals 
      --het \                       
      --out clear

# Remove het fail 
Rscript remove.het.fail.R

# Duplicate snp    
plink --bfile Target_QC \
      --list-duplicate-vars ids-only suppress-first \
      --out dupl
      
plink --bfile Target_QC \
      --exclude dupl.dupvar \
      --make-bed \
      --out target
      
# Check for sex & Impute sex: individuals in which there is a mismatch between biological and reported sex are typically removed
plink --bfile target \
      --extract Target_QC.prune.in \
      --keep clear.het.valid.sample \
      --check-sex \
      --out clear

# Impute sex 
Rscript impute.sex.R

# Relatedness: individuals that have a first or second degree relative in the sample (π>0.125) are removed
plink --bfile target \
      --keep sex.valid \
      --rel-cutoff 0.125 \
      --make-bed \
      --out clear
      
plink --bfile target \
      --make-bed \
      --keep clear.rel.id \
      --set-hh-missing \
      --out Final_Target_QC


#####################################
#                                   #
#             IMPUTATION            #
#                                   #
#####################################

# Download and unzipping the reference genome file
wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

# Converting the PLINK binary files to VCF format
plink --bfile Final_Target_QC \
      --recode vcf \
      --out Target_imputation

# Sorting and compressing the VCF file
bcftools sort Target_imputation.vcf \
        -Oz \
        -o sorted.vcf
bgzip sorted.vcf

# Indexing the sorted VCF files
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
bcftools index -f sorted.vcf.gz 
bcftools index -f All_20180423.vcf.1.gz

# Unzipping the sorted VCF file
gunzip sorted.vcf

# Fixing the reference allele in the VCF file
bcftools +fixref sorted.vcf \
        -Oz \
        -o fixref.vcf \
        -- \
        -d \
        -f Homo_sapiens.GRCh37.dna.primary_assembly.fa \
        -i All_20180423.vcf.1.gz

# Sorting and compressing the fixed reference VCF file
bcftools sort fixref.vcf \
        -Oz \
        -o fixref.sorted.vcf
bgzip fixref.sorted.vcf

# Indexing the sorted VCF files
bcftools index fixref.sorted.vcf.gz -f sorted.vcf.gz

# Splitting the VCF file by chromosome
for i in {1..22}
  do
    bcftools view fixref.sorted.vcf.gz \
              --regions ${i} \
              -o VCF.${i}.vcf.gz \
              -Oz
done


#-----------------------------------------------------------------------------------------------


######     Michigan Imputation Server      ######   
# https://imputationserver.sph.umich.edu/index.html#!run/minimac4
## STEPS
# 1 Select: Run -> Genotype Imputation (Minimac4)
# 2 Reference Panel: HRC r1.1 2016 (GRCh37/hg19)
# 3 Input files: VCF 
# 4 Array Build: GRCh37/hg19
# 5 rsq Filter: NULL
# 6 Phasing: Eagle
# 7 Population: EUR
# 8 Mode: Quality Control & Imputation


#-----------------------------------------------------------------------------------------------


for i in {1..22}
  do
    # Unzipping the specific ZIP file for chromosome
    unzip chr_{i}.zip
    cd /chr_${i}
    
    # Converting the VCF file to PLINK binary format using PLINK2
    plink2 --vcf chr${i}.dose.vcf.gz dosage=HDS \
          --make-bed \
          --out chr
    
    # Applying additional quality control (QC) filters using PLINK on the converted PLINK binary files
    plink --bfile chr \
          --maf 0.01 \
          --hwe 1e-6 \
          --allow-no-sex \
          --write-snplist \
          --make-bed \
          --out Chr
done

# Create: recodedids.txt, recodedids.txt, recodepheno.txt
Rscript information.R

mkdir chr
for i in {1..22}
  do    
    cd chr_${i}
    # Updating sample information based on recoded IDs 
    plink --bfile /Chr \
          --update-ids recodedids.txt \
          --make-bed \
          --out data.id
          
    # Updating sample sex information based on recoded sex IDs 
    plink --bfile data.id \
          --update-sex recodedsex.txt \
          --make-bed \
          --out data.sex
    
    # Adding phenotype information based on recoded phenotype 
    plink --bfile data.sex \
          --pheno recodepheno.txt \
          --make-bed \
          --out data.pheno
    
    # Copying the generated files for each chromosome to the /chr directory
    cp data.pheno.bed /chr/chr${i}.bed
    cp data.pheno.bim /chr/chr${i}.bim
    cp data.pheno.fam /chr/chr${i}.fam
done

cd chr
# Creating a merge list file containing the paths for each chromosome file
for i in {1..22}
  do
    echo chr${i} >> merge.txt
done

# Merging the chromosome files based on merge.txt
plink --merge-list merge.txt \
      --make-bed \
      --out merge
      
#------------------------------------------------------
# How merge the two Target dataset
#-------------------------------------------------------
plink --bfile $path_dataset1/merge 
      --bmerge $path_dataset1/merge
      --make-bed 
      --out merge 
#------------------------------------------------------



#####################################
#                                   #
#                PCA                #
#                                   #
#####################################

# Downloading required files 
wget https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst
wget https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst
wget https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam
      
# Decompressing the downloaded files using PLINK2
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar

# Looping through each chromosome to process the data
for i in {1..22}
  do
    plink2 --pfile all_phase3 vzs \
           --memory 7000 \
           --chr ${i} \
           --snps-only just-acgt \
           --max-alleles 2 \
           --rm-dup exclude-mismatch \
           --set-missing-var-ids '@_#_$1_$2' \
           --new-id-max-allele-len 662 \
           --make-bed \
           --out all.phase3.chr${i}
done

# Creating a merge list file containing the paths for each chromosome file
for i in {1..22}
  do
    echo all.phase3.chr$i >> merge.1000genome.txt
done

# Merging the chromosome files using PLINK based on the merge list file
plink --merge-list merge.1000genome.txt \
      --make-bed \
      --out all.phase3

# Performing filters on the merged file 
plink --bfile merge \
      --snps-only just-acgt \
      --biallelic-only \
      --maf 0.05 \
      --mind 0.1 \
      --hwe 1e-6 \
      --geno 0.05 \
      --make-bed \
      --out merge.db 

plink --bfile merge.db \
      --mind 0.1 \
      --make-bed \
      --out tot.merge.db 
      
# Recovering variant IDs from the merged file using PLINK2
plink2 --bfile tot.merge.db \
       --recover-var-ids all.phase3.bim force partial \
       --make-bed \
       --out Target
       
# Filter reference and study data for non A-T or G-C SNPs
awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") {print $2}' Target.bim > Target.ac_gt_snps
awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") {print $2}' all.phase3.bim > all.phase3.ac_gt_snps

# Create new data sets by excluding the markers identified in the previous step
plink --bfile all.phase3 \
      --exclude all.phase3.ac_gt_snps \
      --make-bed \
      --out all.phase3.no_ac_gt_snps
      
plink --bfile Target \
      --exclude Target.ac_gt_snps \
      --make-bed \
      --out Target.no_ac_gt_snps
      
# Check and correct chromosome mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} ($2 in a && a[$2] != $1) {print a[$2],$2}' Target.no_ac_gt_snps.bim Target.no_ac_gt_snps.bim | sed -n '/^[XY]/!p' > all.phase3.toUpdateChr
plink --bfile all.phase3.no_ac_gt_snps \
      --update-chr all.phase3.toUpdateChr 1 2 \
      --make-bed \
      --out all.phase3.updateChr
      
# Position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} ($2 in a && a[$2] != $4) {print a[$2],$2}' Target.no_ac_gt_snps.bim all.phase3.no_ac_gt_snps.bim > all.phase3.toUpdatePos

# Possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' Target.no_ac_gt_snps.bim all.phase3.no_ac_gt_snps.bim > all.phase3.toFlip

# Upate positions and flip alleles
plink --bfile all.phase3.updateChr \
      --update-map all.phase3.toUpdatePos 1 2 \
      --flip all.phase3.toFlip \
      --make-bed \
      --out all.phase3.flipped
      
# Remove mismatches
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' Target.no_ac_gt_snps.bim all.phase3.autosome.flipped.bim > all.phase3.autosome.mismatch

plink --bfile all.phase3.flipped \
      --exclude all.phase3.mismatch \
      --make-bed \
      --out all.phase3.clean
      
# Pruned
plink --bfile Target.no_ac_gt_snps \
      --indep-pairwise 50 5 0.2 \
      --out Target.no_ac_gt_snps
      
plink --bfile Target.no_ac_gt_snps \
      --extract Target.no_ac_gt_snps.prune.in \
      --make-bed \
      --out Target.pruned
      
plink --bfile all.phase3.clean \
      --extract Target.no_ac_gt_snps.prune.in \
      --make-bed \
      --out all.phase3.pruned
      
# Merge study genotypes and reference data
plink --bfile Target.pruned \
      --bmerge all.phase3.pruned.bed all.phase3.pruned.bim all.phase3.pruned.fam \
      --make-bed \
      --out Target.merge.all.phase3
# PCA
plink --bfile Target.merge.all.phase3 \
      --pca 5 \
      --out data.pca

# PLOT PCA and show the two Target Dataset that overlap
Rscript PLOT_PCA_per_db.R

# Plot PCA of merged dataset
Rscript Plot_PCA.R

# PCA barplot
Rscript PCA_BARPLOT.R

# Remove the people non below to EUR population
plink --bfile Target \
      --keep /newdataset.txt \
      --make-bed \
      --out final.data



#####################################
#                                   #
#             PLINK PRS             #
#                                   #
#####################################

# Add BETA column to base dataset
Rscript clumping.data.R

plink --bfile final.data \
      --clump-p1 1 \
      --clump-r2 0.1 \
      --clump-kb 250 \
      --clump BaseDataset \
      --clump-snp-field SNP \
      --clump-field P \
      --out clumping
      
awk 'NR!=1{print $3}' clumping.clumped >  clumping.valid.snp

# Generate PRS
awk '{print $2,$12}'  BaseDataset > SNP.pvalue
echo "0.001 0 0.001" > range_list 
echo "0.01 0 0.01" >> range_list
echo "0.02 0 0.02" >> range_list
echo "0.03 0 0.03" >> range_list
echo "0.04 0 0.04" >> range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

plink --bfile final.data \
      --score BaseDataset 2 4 15 header \
      --q-score-range range_list SNP.pvalue \
      --extract clumping.valid.snp \
      --out clumping

# Clumping and print the best result
Rscript clumping.fit.R
cat result.txt


#####################################
#                                   #
#             PRSice-2              #
#                                   #
#####################################

#######      PRSice      #######
Rscript PRSice.R

Rscript $path_PRSice/PRSice.R \
        --prsice $path_PRSice/PRSice_linux \
        --base LifetimeMDD \
        --target final.data \
        --binary-target T \
        --pheno phenotype.phen \
        --cov file.covariate \
        --out prs
        
# Grafici
Rscript   Plink.plot.R
Rscript   PRSice.plot.R
Rscript   plink.barplot.R 
