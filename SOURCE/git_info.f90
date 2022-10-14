module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "IntCholesky_AC0block"
    character(len=*), parameter :: GIT_REVISION_HASH = "f722ed4a64d45d4d2fc82697e3eb8c504ea34a28"
    character(len=*), parameter :: GIT_REVISION = "f722ed4"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Mon Oct 10 19:27:48 2022 +0200"
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
