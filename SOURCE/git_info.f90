module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "rdmcorr"
    character(len=*), parameter :: GIT_REVISION_HASH = "fd4b704e2cc4e892c2401e401a568b5c68071cb9"
    character(len=*), parameter :: GIT_REVISION = "fd4b704"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Fri Mar 24 20:46:29 2023 +0100"
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
