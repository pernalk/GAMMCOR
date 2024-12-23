module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "dSRS"
    character(len=*), parameter :: GIT_REVISION_HASH = "8235967b202fa4a8727e7591b67385230cfce9f1"
    character(len=*), parameter :: GIT_REVISION = "8235967"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Sun Dec 22 00:30:29 2024 +0100"
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
