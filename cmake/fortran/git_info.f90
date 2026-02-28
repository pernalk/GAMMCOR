module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "@GIT_HEAD_BRANCH@"
    character(len=*), parameter :: GIT_REVISION_HASH = "@GIT_REVISION_HASH@"
    character(len=*), parameter :: GIT_REVISION = "@GIT_REVISION@"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "@GIT_COMMITTER_DATE@"
    character(len=*), parameter :: GIT_LATEST_TAG_LONG = "@GIT_LATEST_TAG_LONG@"

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
