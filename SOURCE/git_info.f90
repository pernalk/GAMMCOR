module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "IntCholesky"
    character(len=*), parameter :: GIT_REVISION_HASH = "77b93c8b285da236f3b4286c70b1f79414f2f042"
    character(len=*), parameter :: GIT_REVISION = "77b93c8"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Thu Apr 13 16:35:01 2023 +0200"
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
