module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "dSRS"
    character(len=*), parameter :: GIT_REVISION_HASH = "2a78efc386fca1f20f702605add005c3605142aa"
    character(len=*), parameter :: GIT_REVISION = "2a78efc"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Mon Jan 20 17:41:14 2025 +0100"
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
