message(NOTICE "")
message(NOTICE "🏷️   Git info")

# include CPM.cmake module
include(GetCPM)

# add GitInfo.cmake
CPMAddPackage(
    NAME GitInfo.cmake
    GITHUB_REPOSITORY cppmf/GitInfo.cmake
    VERSION 1.0.0
)

# print git info
GitInfo(${CMAKE_CURRENT_SOURCE_DIR})
message("Compiling on branch: .......... ${GIT_HEAD_BRANCH}")
message("Revision (short): ............. ${GIT_REVISION}")
message("Revision (full): .............. ${GIT_REVISION_HASH}")
message("Latest commit date: ........... ${GIT_COMMITTER_DATE}")
message("Lateset tag: .................. ${GIT_LATEST_TAG_LONG}")

# write git info to git_info.f90
configure_file(cmake/fortran/git_info.f90 ${CMAKE_SOURCE_DIR}/SOURCE/git_info.f90)

message(NOTICE "")