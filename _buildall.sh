#!/bin/bash                                                                                                                                   

#cd $(mktemp -d)  #start in a clean directory
#git clone git@github.com:Bioconductor/BiocWorkshops.git
#cd BiocWorkshops
#chmod 755 *.sh


git clone -b individual_builds \
    git@github.com:Bioconductor/BiocWorkshops.git \
    individual_builds

cd individual_builds
git rm -r *.err *.md
cd ..

for file in `ls *.Rmd | grep '^[0-9A-Z]'`; do
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

## copy only .err files containing an error
for file in `find . -type f -name "*.err"`; do
   echo $file
   if grep -q Error $file
   then
      echo $file contains an error
      cp $file individual_builds
   else
      echo OK
   fi
done


cp *.knit.md individual_builds
cd individual_builds
git add --all *.err *.md
git commit -am "individual builds to individual_builds"
git push -q origin master
cd ..

./_build.sh

git clone -b gh-pages \
    git@github.com:Bioconductor/BiocWorkshops.git \
    book-output
cd book-output

git rm -rf *

cp -r ../docs/* ./
git add --all *
git commit -m "Update the book" || true
git push -q origin gh-pages
cd ..
