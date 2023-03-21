module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "IntCholesky"
    character(len=*), parameter :: GIT_REVISION_HASH = "fa5fc92caaa98e37e21e0454b4a69fa8bae80db4"
    character(len=*), parameter :: GIT_REVISION = "fa5fc92"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Fri Feb 17 14:29:13 2023 +0100"
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
