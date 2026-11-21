
module sys_rootdir

implicit none
public

  ! system root directory for ExoRT
  ! Machine-specific by design: this is the ONE file in this bundle that is
  ! expected to differ from source/src.main/ (see the sync rule in
  ! populate3Dmodels.py / CLAUDE.md). Edit the active line for your machine
  ! before building. The path must point at an ExoRT checkout carrying the
  ! v2 data layout (data/kdist/<gas>/, data/stellar/, data/cloud/,
  ! data/aerosol/).

  ! Machine: Summit
  !character(len=256), parameter :: exort_rootdir = '/projects/wolfet/models/ExoRT/'

  ! Machine: Hyak
  !character(len=256), parameter :: exort_rootdir = '/suppscr/vsm/gscratch/wolfet/ExoRT/'

  ! Machine: discover
  ! NOTE: must be a refactor-branch checkout (v2 data layout). The production
  ! /discover/nobackup/etwolf/models/ExoRT checkout is main (v1, data/solar/
  ! layout) and fails the k-file grid check. The first ExoCAM run
  ! (2026-11-21) used /discover/nobackup/etwolf/tmp/exort-bench/ExoRT-refactor/.
  character(len=256), parameter :: exort_rootdir = '/discover/nobackup/etwolf/models/ExoRT/'

end module sys_rootdir
