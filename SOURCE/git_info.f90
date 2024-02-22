module git_info

    implicit none
    character(len=*), parameter :: GIT_HEAD_BRANCH = "dSRS"
    character(len=*), parameter :: GIT_REVISION_HASH = "17e4a2784e9423bd16d464cc63986e146b41d436"
    character(len=*), parameter :: GIT_REVISION = "17e4a27"
    character(len=*), parameter :: GIT_COMMITTER_DATE = "Thu Apr 27 18:23:47 2023 +0200"
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
