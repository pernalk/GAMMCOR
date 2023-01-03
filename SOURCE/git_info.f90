module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "dSRS"
    character(len=*), parameter :: GIT_REVISION_HASH = "863a922a7ec04a5daff86beee26e93a6a493d716"
    character(len=*), parameter :: GIT_REVISION = "863a922"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Wed Dec 28 16:44:47 2022 +0100"
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
