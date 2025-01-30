module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "dSRS"
    character(len=*), parameter :: GIT_REVISION_HASH = "8be98ab74c32592af68da2c83c455e0c7cc45ef0"
    character(len=*), parameter :: GIT_REVISION = "8be98ab"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Wed Jan 22 10:22:58 2025 +0100"
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
