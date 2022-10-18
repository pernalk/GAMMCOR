module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "IntCholesky"
    character(len=*), parameter :: GIT_REVISION_HASH = "3fe3f6ade070f215fe3b1c4ce9d1b13312c64148"
    character(len=*), parameter :: GIT_REVISION = "3fe3f6a"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Tue Oct 18 17:46:32 2022 +0200"
    character(len=*), parameter :: GIT_LATEST_TAG_LONG = "v-4.1"

    contains

    subroutine git_print_info()
        write(*,*) ""
        write(*,*) "GIT INFO"
        write(*,'(a)') "********************************************************************************"
        write(*,*) "Branch:             ", GIT_HEAD_BRANCH
        write(*,*) "Revision (short):   ", GIT_REVISION
        write(*,*) "Revision (full):    ", GIT_REVISION_HASH
        write(*,*) "Latest commit date: ", GIT_COMMITTER_DATE
        write(*,*) "Latest tag:         ", GIT_LATEST_TAG_LONG
        write(*,*) ""
    end subroutine

end module git_info
