module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "aks-openmp"
    character(len=*), parameter :: GIT_REVISION_HASH = "acaa15ddf88ea9b906eefd818402c853d58e50d4"
    character(len=*), parameter :: GIT_REVISION = "acaa15d"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Sun Feb 27 18:08:42 2022 +0100"
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
