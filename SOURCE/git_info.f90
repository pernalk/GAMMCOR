module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "cmake-int-chol"
    character(len=*), parameter :: GIT_REVISION_HASH = "7175d38b90442a488dd29bea3f3b7127273ce6ec"
    character(len=*), parameter :: GIT_REVISION = "7175d38"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Sun Feb 6 02:17:08 2022 +0100"
    character(len=*), parameter :: GIT_LATEST_TAG_LONG = "v-4-🐀"

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
