
module sys_rootdir

implicit none
public

  ! system root directory for ExoRT
  ! KNOWN ISSUE (April, 20206): rootdir is hardcoded path.
  ! edit this line manually before building.
  !
  ! Runtime variable (not a parameter) so the library entry point
  ! exort_init(data_root, ...) can override it without a rebuild.
  ! The standalone executables use the compiled default below.

  ! Machine: Summit
  !character(len=256) :: exort_rootdir = '/projects/wolfet/models/ExoRT/'

  ! Machine: Hyak
  !character(len=256) :: exort_rootdir = '/suppscr/vsm/gscratch/wolfet/ExoRT/'

  ! Machine: discover
  !character(len=256) :: exort_rootdir = '/discover/nobackup/etwolf/models/ExoRT/'

  character(len=256) :: exort_rootdir = '/Users/wolfe/models/ExoRT/'

end module sys_rootdir
