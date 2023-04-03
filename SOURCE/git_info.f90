module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "IntCholesky"
    character(len=*), parameter :: GIT_REVISION_HASH = "3873c1f17f00ed54ec5a879f598782f8aef42919"
    character(len=*), parameter :: GIT_REVISION = "3873c1f"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Wed Mar 29 18:57:54 2023 +0200"
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
