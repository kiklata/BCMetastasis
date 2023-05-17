for file in /home/zhepan/BCLM/down/*.fastq.gz;
do
  if [ -f "$file" ]; then
    samples=${file:39:4}
    newdir="/home/zhepan/BCLM/fastq/NGS$samples"
    mkdir -p "$newdir"
    newfile="$newdir/$(basename "$file")"
    printf "Moving %s to %s\n" "$file" "$newfile"
    mv -i "$file" "$newfile"
  fi
done
