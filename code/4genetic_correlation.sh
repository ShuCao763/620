#!/bin/bash
# 设置路径
LDSC_PATH="/Users/caoshu/ldsc-2.0.1"
RESULT_DIR="/Users/caoshu/WISC/620/final_project/res_for_genecorr"
MUNGED_DIR="/Users/caoshu/WISC/620/final_project/munged_data"
REF_LD_CHR="/Users/caoshu/WISC/620/ldsc_inputs/for_h2/eur_w_ld_chr/"
W_LD_CHR="/Users/caoshu/WISC/620/ldsc_inputs/for_enrichment/weights/weights.hm3_noMHC."

# 定义目标性状和比较性状
target_trait="Chronotype"
compare_traits=("ADHD" "EA" "Baldness" "SCZ" "Napping")

# 创建输出汇总表头（如果不存在）
summary_file="$RESULT_DIR/genetic_correlation_summary.tsv"
if [ ! -f "$summary_file" ]; then
  echo -e "Trait1\tTrait2\trg\trg_se\tz\tp" > "$summary_file"
fi

# 遍历并计算与 chronotype 的 correlation
for trait in "${compare_traits[@]}"; do
  out_prefix="$RESULT_DIR/genetic_correlation_${target_trait}_${trait}"
  echo "Running LDSC rg: $target_trait vs $trait"

  python $LDSC_PATH/ldsc.py \
    --rg "$MUNGED_DIR/munged_${target_trait}.sumstats.gz,$MUNGED_DIR/munged_${trait}.sumstats.gz" \
    --ref-ld-chr "$REF_LD_CHR" \
    --w-ld-chr "$W_LD_CHR" \
    --out "$out_prefix"

  logfile="${out_prefix}.log"
  summary_line=$(awk '/^Summary of Genetic Correlation Results/ {getline; getline; print}' "$logfile")

  # 提取字段
  p1=$(echo "$summary_line" | awk '{print $1}')
  p2=$(echo "$summary_line" | awk '{print $2}')
  rg=$(echo "$summary_line" | awk '{print $3}')
  se=$(echo "$summary_line" | awk '{print $4}')
  z=$(echo "$summary_line" | awk '{print $5}')
  p=$(echo "$summary_line" | awk '{print $6}')

  trait1=$(basename "$p1" .sumstats.gz | sed 's/^munged_//')
  trait2=$(basename "$p2" .sumstats.gz | sed 's/^munged_//')

  echo -e "$trait1\t$trait2\t$rg\t$se\t$z\t$p" >> "$summary_file"
done

#最后直接 输入： bash 4genetic_correlation.sh
