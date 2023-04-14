module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "workshop2023"
    character(len=*), parameter :: GIT_REVISION_HASH = "e6323b4b44c8efb791e030ba8cd2ca2c10815dab"
    character(len=*), parameter :: GIT_REVISION = "e6323b4"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Fri Apr 14 13:36:51 2023 +0200"
    character(len=*), parameter :: GIT_LATEST_TAG_LONG = ""

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
