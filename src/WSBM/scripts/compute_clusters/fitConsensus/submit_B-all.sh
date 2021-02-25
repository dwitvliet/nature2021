folder=$(echo $0 | cut -d '_' -f2 | cut -d '.' -f1)

for i in 0 1 4 6
do
  echo submit $i
  sbatch $folder/submit_$i.sh
done

