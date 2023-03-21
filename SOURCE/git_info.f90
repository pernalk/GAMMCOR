module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "rdmcorr"
    character(len=*), parameter :: GIT_REVISION_HASH = "373453a64b4bccd32d3898b5bc270df4bcfa65c6"
    character(len=*), parameter :: GIT_REVISION = "373453a"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Fri Mar 10 16:51:19 2023 +0100"
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
