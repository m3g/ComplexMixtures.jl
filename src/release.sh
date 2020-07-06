#!/bin/bash
####################################################################################################

# Software name:

package=MDDF

# GIT URL:

giturl=https://github.com/m3g/MDDF

# Name of file containing version number

versionfile=./version.jl

####################################################################################################

#git log --pretty=oneline 16.323...16.330 | awk '{$1=""; print "-"$0}'

year=`date +%y`
day=`date +%j`
baseversion="${year:0:1}${year:1:1}.$day"

taglist=`git tag --list "$baseversion"`

# While the current tag list contains the version, increase version number
version="$baseversion"
i=0
while [[ "$taglist" == *"$version"* ]] ; do
  i=`expr $i + 1`
  version="$baseversion.$i"
done

echo " Tagging with version: $version" 

cat $versionfile | sed -e "s/Version.*/Version\ $version\"/" > versiontmp
\mv -f versiontmp $versionfile

git add -A 
git commit -m "Changed version file to $version"
git tag -a $version -m "Release $version"
git push origin master tag $version
 
echo "----------------------"
echo "CHANGE LOG:"
echo "----------------------"
range=`git tag | tail -n 2 | xargs | sed 's! !...!'`
git log --pretty=oneline $range | awk '{$1=""; print "-"$0}'
echo "----------------------"

echo " Done. " 

