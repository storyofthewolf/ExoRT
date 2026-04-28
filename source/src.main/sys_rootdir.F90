
module sys_rootdir

implicit none
public

  ! system root directory for ExoRT
  ! KNOWN ISSUE (April, 20206): rootdir is hardcoded path.
  ! edit this line manually before building.

  ! Machine: Summit
  !character(len=256), parameter :: exort_rootdir = '/projects/wolfet/models/ExoRT/'  

  ! Machine: Hyak
  !character(len=256), parameter :: exort_rootdir = '/suppscr/vsm/gscratch/wolfet/ExoRT/'  

  ! Machine: discover
  !character(len=256), parameter :: exort_rootdir = '/discover/nobackup/etwolf/models/ExoRT/'

  character(len=256), parameter :: exort_rootdir = '/Users/wolfe/models/ExoRT/'

end module sys_rootdir
