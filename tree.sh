#!/bin/bash
#####################################################
# Unix/Linux Folder Tree                            #
# Version: 2.8                                      #
# File: ~/apps/tree/tree.sh                         #
#                                                   #
# Displays structure of folder hierarchy            #
# ------------------------------------------------- #
# This tiny script uses "ls", "grep", and "sed"     #
# in a single command to show the nesting of        #
# subfolders.                                       #
#                                                   #
# Setup:                                            #
#    $ cd ~/apps/tree                               #
#    $ chmod +x tree.sh                             #
#    $ ln -s tree.sh /usr/local/bin/tree            #
#    $ which tree                                   #
#                                                   #
# Usage:                                            #
#    $ tree [FOLDER]                                #
#                                                   #
# Examples:                                         #
#    $ tree                                         #
#    $ tree /etc/apache2                            #
#    $ tree ..                                      #
#                                                   #
# WTFPL ~ https://centerkey.com/tree ~ Dem Pilafian #
#####################################################

echo
test -z "$1" || cd "$1"  #if parameter exists, use as base folder
pwd
ls -R | grep "^[.]/" | sed -e "s/:$//" -e "s/[^\/]*\//--/g" -e "s/^/   |/" 
# grep:    select folders (filter out files)
# 1st sed: remove trailing colon
# 2nd sed: replace higher level folder names with dashes
# 3rd sed: indent graph and add leading vertical bar
topFolders=$(ls -F -1 | grep "/" | wc -l)
test $topFolders -ne 0 || echo "   --> no subfolders"
echo