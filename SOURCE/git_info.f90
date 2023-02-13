module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "IntOnTheFly"
    character(len=*), parameter :: GIT_REVISION_HASH = "0c8f28b7d73c7331da2b2a1e46e00768e41c2be3"
    character(len=*), parameter :: GIT_REVISION = "0c8f28b"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Fri Nov 25 12:12:24 2022 +0100"
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
