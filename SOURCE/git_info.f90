module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "IntOnTheFly"
    character(len=*), parameter :: GIT_REVISION_HASH = "e00309315a9dfa3715206aaacde8346d51543ea9"
    character(len=*), parameter :: GIT_REVISION = "e003093"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Tue Jul 25 10:38:50 2023 +0200"
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
