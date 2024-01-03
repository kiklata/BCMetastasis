path=/home/shpc_100839/PaperCD8/data/Tex/compass/liu
outpath=/home/shpc_100839/PaperCD8/data/Tex/compass/res

ls $path/*.tsv|cut -d"/" -f 9|cut -d"." -f 1 | sort -u |
while read id;do 
  echo ${id}
  mkdir $outpath/${id} 
  compass --data $path/${id}.tsv --num-processes 20 --species homo_sapiens --output-dir $outpath/${id} 
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done
