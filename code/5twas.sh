#!/bin/bash
#用自己的电脑，不用server

#1 clone文件和下载必要文件
git clone https://github.com/hakyimlab/MetaXcan
https://uchicago.app.box.com/s/us7qhue3juubq66tktpogeansahxszg9（到这下载数据）

#2.1 调环境1(输入下面两行，如果不可以，用2.2中的方式调环境)
conda env create -f /path/to/this/repo/software/conda_env.yaml
conda activate imlabtools  

#2.2 调环境2（必须按下面的一行行输入）
CONDA_SUBDIR=osx-64 conda create -n imlabtools_env python=3.7
conda activate imlabtools_env
CONDA_SUBDIR=osx-64 conda install numpy=1.18.1 pandas=0.25.3 scipy=1.4.1 statsmodels=0.11.1 pyarrow=0.11.0 -c conda-forge
CONDA_SUBDIR=osx-64 conda install cyvcf2=0.20.0 bgen_reader=3.0.2 -c bioconda -c conda-forge

conda activate imlabtools_env
#3（进入conda环境后，terminal光标会显示‘imlabtools_env’，然后直接就可以跑了）
python /Users/caoshu/WISC/620/MetaXcan/software/SPrediXcan.py \
--model_db_path /Users/caoshu/WISC/620/MetaXcan/example_data/en_Brain_Hypothalamus.db \
--covariance /Users/caoshu/WISC/620/MetaXcan/example_data/en_Brain_Hypothalamus.txt.gz \
--gwas_file /Users/caoshu/WISC/620/final_project/data/Chronotype.txt.gz \
--gwas_N 449734 \
--snp_column SNP \
--effect_allele_column ALLELE1 \
--non_effect_allele_column ALLELE0 \
--beta_column BETA \
--pvalue_column P_BOLT_LMM \
--output_file /Users/caoshu/WISC/620/final_project/result/twas3.csv

conda deactivate

#测试部分（忽略）
python /Users/caoshu/WISC/620/MetaXcan/software/SPrediXcan.py \
--model_db_path /Users/caoshu/WISC/620/MetaXcan/example_data/DGN-WB_0.5.db \
--covariance /Users/caoshu/WISC/620/MetaXcan/example_data/covariance.DGN-WB_0.5.txt.gz \
--gwas_folder /Users/caoshu/WISC/620/MetaXcan/example_data/GWAS \
--gwas_file_pattern ".*gz" \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column BETA \
--pvalue_column P \
--output_file /Users/caoshu/WISC/620/final_project/result/twas_test.csv

#整体数据(按P值)
twas <- fread("/Users/caoshu/WISC/620/final_project/result/twas3.csv")
twas_sorted <- twas[order(twas$pvalue), ]
fwrite(twas_sorted, "/Users/caoshu/WISC/620/final_project/result/twas_brain_hypoth.csv", row.names = FALSE)
