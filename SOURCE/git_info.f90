module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "IntOnTheFly"
    character(len=*), parameter :: GIT_REVISION_HASH = "200dce0b5dbff4f1ac5fac9f76aee9e60d3d896b"
    character(len=*), parameter :: GIT_REVISION = "200dce0"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Tue Sep 26 14:40:29 2023 +0200"
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
