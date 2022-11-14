module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "IntOnTheFly"
    character(len=*), parameter :: GIT_REVISION_HASH = "ab0a6ff09d229c1f0a91fbd80ae3dda44b6e9f6f"
    character(len=*), parameter :: GIT_REVISION = "ab0a6ff"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Tue Nov 8 19:09:09 2022 +0100"
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
