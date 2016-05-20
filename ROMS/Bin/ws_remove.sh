#!/bin/csh -f
#
# svn $Id: ws_remove.sh 429 2009-12-20 17:30:26Z arango $
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
#    ws_remove.sh [options]                                          :::
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

set log = 0
set verb = 0
set all = -1
set c_dirs = ""
set special_files = ""
set tab_list = ""
set s_count = 0
set t_count = 0

# Set spaces and tabs log file names

set spaces = "spaces_removed.log"
set tabs = "tabs_found.log"

while ( ($#argv) > 0 )
  switch ($1)

   case "-v":
      shift
      set verb = 1
    breaksw

    case "--log":
      shift
      set log = 1

      # Remove the log files
      /bin/rm -f $spaces $tabs
    breaksw

     case "--all":
      shift

      if ( $all == 0 ) then
        echo "-f and --all cannot be used together"
        exit 1
      endif

      set c_dirs = "Compilers Master ROMS User"
      set special_files = "makefile Waves/SWAN/Src/Module.mk Waves/SWAN/Src/waves_coupler.F"
      set all = 1
    breaksw

    case "-f":
      shift

      if ( $all == 1 ) then
        echo "-f and --all cannot be used together"
        exit 1
      endif

      @ c = 0
      set all = 0
      while ( `echo $1 | grep -c '^-'` == 0 && $1 != "" )
        set special_files = "${special_files}${1} "
        @ c++
        shift
      end
      if ( $c == 0 ) then
        echo "-f must be followed by filename(s)"
        exit 1
      endif
    breaksw

   case "-*":
      set cmd = `basename $0`
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
    breaksw
  endsw
end

if ( $all == -1 ) then
  set cmd = `basename $0`
  echo ""
  echo "${cmd}: You must use at least the -f or --all option"
  echo ""
  exit 1
endif

if ( $all == 1 ) then
  /bin/echo -e "\nGenerating file list ...\n"
else
  /bin/echo -e "\nThe following file(s) will be processed:\n"
  /bin/echo -e "${special_files}\n"
endif

# Set space and tab counters (number of files not number of occurances)

@ s=0
@ t=0

# The "! -path '*/.svn/*'" is there to keep it from messing with
# files in the .svn directories. There is no way to redirect only
# stderr with csh

if ( $all == 1 ) then
  foreach FILE ( `find ${c_dirs} ! -path '*/.svn/*' -type f -print | sort` )
    set s_count = 0
    set t_count = 0
    # Double check that we're not changing a file in a .svn folder
    if ( `/bin/echo $FILE | grep -vc '.svn/'` ) then
      # Grep for trailing white space including tabs
#      set s_count = `grep -cP '[ \t]+$' $FILE`
      set s_count = `grep -c '[ \t]+$' $FILE`
      if ( $s_count >= 1 ) then
        # Increment spaces file counter
        @ s++

        # Add FILE to spaces log if --log set
        if ( $log == 1 ) then
          /bin/echo -en "${FILE} ... " >> $spaces
        endif

        # Also output to screen if -v is set
        if ( $verb == 1 ) then
          /bin/echo -en "$FILE ... "
          # Actual replacement command
          perl -ei 's|[ \t]*$||g' $FILE
#          sed -i 's|[ \t]*$||g' $FILE
          /bin/echo -en "${s_count} replacement(s) made\n"
        else
          # Actual replacement command
          perl -ei 's|[ \t]*$||g' $FILE
#          sed -i 's|[ \t]*$||g' $FILE
        endif

        # Report in spaces log that FILE has been processed if --log set
        if ( $log == 1 ) then
          /bin/echo -en "${s_count} replacement(s) made\n" >> $spaces
        endif
      endif

      # make list of files with tabs to edit by hand
      set t_count = `grep -c "\t" $FILE`
#      set t_count = `grep -cP "\t" $FILE`
      if ( $t_count >= 1 ) then
        # Increment tabs file counter and add FILE to tabs log
        @ t++

        # Add FILE to tabs log if --log set; otherwise update $tab_list
        if ( $log == 1 ) then
          /bin/echo "${FILE} ... ${t_count} tab(s) found" >> $tabs
        endif
        set tab_list = "${tab_list}${FILE} ... ${t_count} tab(s) found\n"
      endif

    else
      /bin/echo -e "\nThere is a .svn in the path: $FILE skipped\n"
    endif
  end
endif

# Now look in the special_files

foreach FILE ( $special_files )
  set s_count = 0
  set t_count = 0
  # Grep for trailing white space including tabs
  set s_count = `grep -c '[ \t]+$' $FILE`
#  set s_count = `grep -cP '[ \t]+$' $FILE`
  if ( $s_count >= 1 ) then
    # Increment spaces file counter
    @ s++

    # Add FILE to spaces log if --log set
    if ( $log == 1 ) then
      /bin/echo -en "${FILE} ... " >> $spaces
    endif

    # Also output to screen if -v is set
    if ( $verb == 1 ) then
      /bin/echo -en "${FILE} ... "
      # Actual replacement command
      perl -ei 's|[ \t]*$||g' $FILE
#      sed -i 's|[ \t]*$||g' $FILE
      /bin/echo -en "${s_count} replacement(s) made\n"
    else
      # Actual replacement command
      perl -ei 's|[ \t]*$||g' $FILE
#      sed -i 's|[ \t]*$||g' $FILE
    endif

    # Report in spaces log that FILE has been processed if --log set
    if ( $log == 1 ) then
      /bin/echo -en "${s_count} replacement(s) made\n" >> $spaces
    endif
  endif

  # Finish tabs file list. We need to leave the tabs in the
  # makefile though.
  set t_count = `grep -c '\t' $FILE`
#  set t_count = `grep -cP '\t' $FILE`
  if ( $t_count >= 1 && $FILE != "makefile" ) then
    # Increment tabs file counter and add FILE to tabs log
    @ t++

    # Add FILE to tabs log if --log set and update $tab_list
    if ( $log == 1 ) then
      /bin/echo "${FILE} ... ${t_count} tab(s) found" >> $tabs
    endif
    set tab_list = "${tab_list}${FILE} ... ${t_count} tab(s) found\n"
  endif
end

/bin/echo -e "\n# of files with trailing spaces: $s"
/bin/echo "# of files with tabs:            $t"

if ( $t >= 1 ) then
  /bin/echo -e "\n********************************************************"
  /bin/echo -e "********************************************************"
  /bin/echo -e "**                                                    **"
  /bin/echo -e "** The files listed below should be checked by hand   **"
  /bin/echo -e "** to remove the remaining tabs.                      **"
  /bin/echo -e "**                                                    **"
  /bin/echo -e "********************************************************"
  /bin/echo -e "********************************************************\n"
  /bin/echo -en "$tab_list"
endif

/bin/echo -e "\nFinished.\n"
