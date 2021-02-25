folder=$(echo $0 | cut -d '_' -f2 | cut -d '.' -f1)

for i in 1
do
  echo "$folder/submit_$i.sh"
  sbatch $folder/submit_$i.sh
  sleep 1
done

