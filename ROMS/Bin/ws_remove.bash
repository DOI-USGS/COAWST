#!/bin/bash
#
# svn $Id: ws_remove.bash 429 2009-12-20 17:30:26Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2016 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS/TOMS White Space Removal Script                                  :::
#                                                                       :::
# Script to remove trailing white space (including tabs) and report all :::
# files that contain tabs (excluding the makefile) so you can remove    :::
# the tabs by hand.                                                     :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ws_remove.bash [options]                                           :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    --all      remove trailing spaces from all ROMS files              :::
#                                                                       :::
#    -f file(s) remove trailing spaces from listed files (wildcards ok) :::
#                                                                       :::
#    --log      output files that are modified to .log files            :::
#                                                                       :::
#    -v         output files that are modified to screen                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Initialize variables

log=0
verbose=0
all=-1
c_dirs=""
special_files=""
tab_list=""
s_count=0
t_count=0

# Set spaces and tabs log file names

spaces="spaces_removed.log"
tabs="tabs_found.log"

while [ $# -gt 0 ]
do
  case "$1" in
    -v )
      shift
      verbose=1
    ;;

    --log )
      shift
      log=1

      # Remove the log files
      /bin/rm -f $spaces $tabs
    ;;

    --all )
      shift

      if [ $all -eq 0 ]; then
        echo "-f and --all cannot be used together"
        exit 1
      fi

      c_dirs="Compilers Master ROMS User"
      special_files="makefile Waves/SWAN/Src/Module.mk Waves/SWAN/Src/waves_coupler.F"
      all=1
    ;;

    -f )
      shift

      if [ $all -eq 1 ]; then
        echo "-f and --all cannot be used together"
        exit 1
      fi

      c=0
      all=0
      #if [ $1 ]; then
      while [[ $1 != '' ]] && [ `echo $1 | grep -c '^-'` -eq 0 ]
      do
        special_files="${special_files}${1} "
        let c=c+1
        shift
      done
      #fi
      if [ $c -eq 0 ]; then
        echo "-f must be followed by filename(s)"
        exit 1
      fi
    ;;

    * )
      cmd=`basename $0`
      echo ""
      echo "${cmd} : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "   --all       remove trailing spaces from all ROMS files"
      echo ""
      echo "   -f file(s)  remove trailing spaces from listed files (wildcards ok)"
      echo ""
      echo "   --log       output files that are modified to .log files"
      echo ""
      echo "   -v          output files that are modified to screen"
      echo ""
      exit 1
    ;;
  esac
done

if [ $all -eq -1 ]; then
  cmd=`basename $0`
  echo ""
  echo "${cmd}: You must use at least the -f or --all option"
  echo ""
  exit 1
fi

if [ $all -eq 1 ]; then
  echo -e "\nGenerating file list ...\n"
else
  echo -e "\nThe following file(s) will be processed:\n"
  echo -e "${special_files}\n"
fi

# Set space and tab counters (number of files not number of occurances)

s=0
t=0

# The "! -path '*/.svn/*'" is there to keep it from messing with
# files in the .svn directories. "2>" redirects stderr so errors
# don't get put in FILE

if [ $all -eq 1 ]; then
  for FILE in `find ${c_dirs} ! -path '*/.svn/*' -type f -print 2> /dev/null | sort`
  do
    s_count=0
    t_count=0
    # Double check that we're not changing a file in a .svn folder
    if [ `echo $FILE | grep -vc '.svn/'` -gt 0 ]; then
      # Grep for trailing white space including tabs
      s_count=`grep -cP "[ \t]+$" $FILE`
      if [ $s_count -ge 1 ]; then
        # Increment spaces file counter
        let s=s+1

        # Add FILE to spaces log if --log set
        if [ $log -eq 1 ]; then
          echo -en "${FILE} ... " >> $spaces
        fi

        # Also output to screen if -v is set
        if [ $verbose -eq 1 ]; then
          echo -en "${FILE} ... "
          # Actual replacement command
          sed -i 's|[ \t]*$||g' $FILE
          echo -en "${s_count} replacement(s) made\n"
        else
          # Actual replacement command
          sed -i 's|[ \t]*$||g' $FILE
        fi

        # Report in spaces log that FILE has been processed if --log set
        if [ $log -eq 1 ]; then
          echo -en "$s_count replacement(s) made\n" >> $spaces
        fi
      fi

      # make list of files with tabs to edit by hand
      t_count=`grep -cP "\t" $FILE`
      if [ $t_count -ge 1 ]; then
        # Increment tabs file counter
        let t=t+1

        # Add FILE to tabs log if --log set and update $tab_list
        if [ $log -eq 1 ]; then
          echo "$FILE ... ${t_count} tab(s) found" >> $tabs
        fi
        tab_list="${tab_list}${FILE} ... ${t_count} tab(s) found\n"
      fi

    else
      echo -e "\nThere is a .svn in the path: $FILE skipped\n"
    fi
  done
fi

# Now look in the special_files

for FILE in $special_files
do
  s_count=0
  t_count=0
  # Grep for trailing white space including tabs
  s_count=`grep -cP "[ \t]+$" $FILE`
  if [ $s_count -ge 1 ]; then
    # Increment spaces file counter
    let s=s+1

    # Add FILE to spaces log if --log set
    if [ $log -eq 1 ]; then
      echo -en "${FILE} ... " >> $spaces
    fi

    # Also output to screen if -v is set
    if [ $verbose -eq 1 ]; then
      echo -en "${FILE} ... "
      # Actual replacement command
      sed -i 's|[ \t]*$||g' $FILE
      echo -en "$s_count replacement(s) made\n"
    else
      # Actual replacement command
      sed -i 's|[ \t]*$||g' $FILE
    fi

    # Report in spaces log that FILE has been processed if --log set
    if [ $log -eq 1 ]; then
      echo -en "$s_count replacement(s) made\n" >> $spaces
    fi
  fi

  # Finish tabs file list. We need to leave the tabs in the
  # makefile though.
  t_count=`grep -cP "\t" $FILE`
  if [ $t_count -ge 1 ] && [[ $FILE != "makefile" ]]; then
    # Increment tabs file counter
    let t=t+1

    # Add FILE to tabs log if --log set; otherwise update $tab_list
    if [ $log -eq 1 ]; then
      echo "$FILE ... ${t_count} tab(s) found" >> $tabs
    fi
    tab_list="${tab_list}${FILE} ... ${t_count} tab(s) found\n"
  fi
done

echo -e "\n# of files with trailing spaces: $s"
echo "# of files with tabs:            $t"

if [ $t -ge 1 ]; then
  echo -e "\n********************************************************"
  echo -e "********************************************************"
  echo -e "**                                                    **"
  echo -e "** The files listed below should be checked by hand   **"
  echo -e "** to remove the remaining tabs.                      **"
  echo -e "**                                                    **"
  echo -e "********************************************************"
  echo -e "********************************************************\n"
  echo -en "$tab_list"
fi

echo -e "\nFinished.\n"
