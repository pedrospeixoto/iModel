#!/bin/bash
# P. Peixoto - Jul 2012
# modified by Luan Santos - 2022

date=` date +%F `
version=` date +%y.%m.%d `

output="imodel.tar.bz2"

#Local host backup directory
dropdir="/home/jeferson/Dropbox/doc/code/imodel"

#-------------------------------------------------------------------------------------------------------
# remote host 1 - ime.usp.br
user_remote_host1="Jeferson-local"
remote_host1="brucutuiv.ime.usp.br"
remote_host1_dir="/var/tmp/jbram"
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# remote host 2 - ybytu
user_remote_host2="jbrambatti"
remote_host2="ybytu.ime.usp.br"
remote_host2_dir="/home/jbrambatti/imodel"
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# Send data to local host
echo "Sync with Dropbox:"
rsync -v -t -u $output  "$dropdir/."
echo "Synchronized with Dropbox"
echo
#-------------------------------------------------------------------------------------------------------
 
#-------------------------------------------------------------------------------------------------------
#remote server ime.usp.br backup sync
echo "Sending to $remote_host1:"
rsync -t -v -z -a --progress $output $user_remote_host1@$remote_host1:$remote_host1_dir
echo "Sent to $user_remote_host1@$remote_host1"
echo
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#remote server ybytu backup sync
echo "Sending to $remote_host2:"
ssh -t $user_remote_host1@$remote_host1 "rsync -t -v -z -a --progress $remote_host1_dir/$output $user_remote_host2@$remote_host2:$remote_host2_dir; rm -rf $output"
echo "Sent to $user_remote_host2@$remote_host2"
echo
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# untar and compile
echo "Untar and compilation at $remote_host2"
ssh -t $user_remote_host1@$remote_host1 "ssh -t $user_remote_host2@$remote_host2 <<EOF
	cd $remote_host2_dir;
	tar -xvf $output;
	make clean;
	make;
	rm -rf $output;
EOF"
echo "Untar and compilation at $remote_host2 done."
#-------------------------------------------------------------------------------------------------------

# remove tar file
rm -rf $output
