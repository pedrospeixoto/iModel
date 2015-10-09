#!/bin/bash

# Creates a backup sync of main model files
# P. Peixoto - Jul 2012

date=` date +%F `
version=` date +%y.%m.%d `
#echo "Today: " $date

source src/tarfiles.sh

output="imodel$version.tar.bz2"

OS=`uname`
echo "System:" $OS
echo

#Edit place to sync relative to system used
# It is set to a dropbox directory
#if [ "$OS" == "Linux" ] ; then    
#    dropdir="/home/pedrosp/Dropbox/Work/imodel"
#else
#    dropdir="/cygdrive/c/Documents\ and\ Settings/Pedro/Dropbox/Doutorado/imodel"
#fi
#echo $dropdir/src
#echo "Sync with Dropbox:"
#rsync -v -t -u $output  "$dropdir/."
#echo "Synchronized with Dropbox"
#echo

#Edit remote server backup sync

#echo "Sending to labmap2.ime.usp.br:"
#rsync -t -v -z -a --progress $output pedrosp@labmap2.ime.usp.br:imodel
#echo "Sent to pedrosp@labmap2.ime.usp.br:imodel"

echo "Sending to ime.usp.br:"
rsync -t -v -z -a --progress $output pedrosp@ime.usp.br:imodel
echo "Sent to pedrosp@ime.usp.br:imodel"
echo

#Remove tar file to a backup directory
mv $output ../../Backups/tmp/.
echo "Tar file saved in backups folder"
echo
echo "-------------------------------------------"
