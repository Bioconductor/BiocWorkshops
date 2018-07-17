#!/bin/bash                                                                                                                                   

i=0
for file in `ls *.Rmd | grep '^[0-9A-Z]'`; do
  i=i++
  log=${file%.Rmd}.log
  err=${file%.Rmd}.err
  tmpdir=${file%.Rmd}_tmp
  mkdir -p $tmpdir
  cd $tmpdir
  ln -f -s ../$file
  nice -n 20 Rscript -e "rmarkdown::render('$file', clean=FALSE, params = list(cache=TRUE))" > $log 2> $err &
  cd ..
done

## Wait for all jobs to finish
for job in `jobs -p`
do
echo $job
    wait $job || let "FAIL+=1"
done

echo $FAIL

mkdir -p individual_md
rm individual_md/*.err
rm individual_md/*.md

## copy only .err files containing an error
for file in `find . -type f -name "*.err"`; do
   echo $file
   if grep -q Error $file
   then
      echo $file contains an error
      cp $file individual_md
   else
      echo OK
   fi
done


cp *.knit.md individual_md
git stage individual_md/
git commit -m "individual builds to individual_md"
git push